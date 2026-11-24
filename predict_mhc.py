#!/usr/bin/env python3
"""
Predict peptide binding to murine MHC molecules.

Peptide length in wt_short_peptide / mt_short_peptide determines the MHC class:

    Class I  (lengths 8–12) -> NetMHCpan 4.2 + MixMHCpred 3.0
                               alleles: both H-2-Kb and H-2-Db
    Class II (length 15)    -> NetMHCIIpan 4.3 + hmmhc + MixMHC2pred
                               allele: H-2-IAb (MixMHC2pred is optional)

Every class-I peptide is evaluated against both Kb and Db alleles.

Results for each allele are written to a separate <filename>.<allele>.tsv file.
Each output contains the original columns plus predictions for that allele,
with columns named using the {wt|mt}_{tool}_{metric} convention. Input files
are never modified.

Peptides with unsupported lengths are skipped. Multiple input files are supported.

Example:
    python predict_mhc.py \
        -i file1.tsv file2.tsv \
        --netmhcpan /path/to/netMHCpan \
        --mixmhcpred /path/to/MixMHCpred \
        --netmhc2pan /path/to/netMHCIIpan \
        --hmmhc /path/to/hmmhc-predict

Requires Python 3 and pandas. External prediction tools must be installed
separately and made available on PATH or supplied using their command-line options.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


PEPTIDE_COLUMNS = {"wt": "wt_short_peptide", "mt": "mt_short_peptide"}
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Peptide lengths mapped to MHC class
CLASS_I_LENGTHS = {8, 9, 10, 11, 12}
CLASS_II_LENGTHS = {15}

#   NetMHCpan:   netMHCpan -listMHC | grep H-2
#   MixMHCpred:  check the allele list in the installation directory (lib/alleles...); names may vary by version
#   NetMHCIIpan: netMHCIIpan -list | grep -i IA
# NOTE: tools use different naming conventions for the same alleles.
#   NetMHCpan/NetMHCIIpan: H-2-Kb, H-2-Db, H-2-IAb
#   MixMHCpred 3.0:        H2-Kb, H2-Db (output columns may use Score_H2-Kb / %Rank_H2-Kb)
#   hmmhc:                 H2-IAb
NET_ALLELES   = ["H-2-Kb", "H-2-Db"]   # NetMHCpan, class I
MIX_ALLELES   = ["H2-Kb", "H2-Db"]     # MixMHCpred, class I
NET2_ALLELES  = ["H-2-IAb"]            # NetMHCIIpan, class II
HMMHC_ALLELES = ["H2-IAb"]             # hmmhc, class II
MIX2_ALLELES  = ["H-2-IAb"]            # MixMHC2pred, class II (verify against Alleles_list_*.txt)

# Write a separate output file for each allele.
# Each output combines columns from the listed (tool, tool-specific allele) pairs.
# Allele names may differ between tools even when they refer to the same allele.
OUTPUT_GROUPS = {
    "H-2-Kb":  [("netmhcpan", "H-2-Kb"), ("mixmhcpred", "H2-Kb")],
    "H-2-Db":  [("netmhcpan", "H-2-Db"), ("mixmhcpred", "H2-Db")],
    "H-2-IAb": [("netmhc2pan", "H-2-IAb"), ("hmmhc", "H2-IAb"), ("mixmhc2pred", "H-2-IAb")],
}

# Leading (non-allele-specific) and auxiliary columns in -xls output
LEAD_COLS = {"pos", "peptide", "id", "target", "n"}
SKIP_COLS = {"core", "icore", "of", "core_rel", "exp_bind", "ave", "nb", "pos", "id", "target"}

OLD_COL_RE = re.compile(r"^(wt|mt)_(netmhcpan|netmhc2pan|mixmhcpred|hmmhc)_")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", nargs="+", required=True,
                   help="One or more input TSV files")
    p.add_argument("--netmhcpan", default="netMHCpan", help="Path to NetMHCpan 4.2 (class I)")
    p.add_argument("--mixmhcpred", default="MixMHCpred", help="Path to MixMHCpred 3.0 (class I)")
    p.add_argument("--netmhc2pan", default="netMHCIIpan", help="Path to NetMHCIIpan 4.3 (class II)")
    p.add_argument("--hmmhc", default="hmmhc-predict", help="Path to hmmhc-predict (class II)")
    p.add_argument("--mixmhc2pred", default=None,
                   help="Path to MixMHC2pred (class II); omitted unless explicitly provided")
    p.add_argument("--mixmhc2pred-context", action="store_true",
                   help="Pass context to MixMHC2pred (default: --no_context)")
    p.add_argument("--no-ba", action="store_true", help="Do not request binding-affinity output from NetMHC* (faster)")
    p.add_argument("--skip-mixmhcpred", action="store_true")
    p.add_argument("--skip-netmhc2pan", action="store_true")
    p.add_argument("--skip-hmmhc", action="store_true")
    p.add_argument("--skip-mixmhc2pred", action="store_true")
    return p.parse_args()


def check_executable(path, name):
    resolved = shutil.which(path) or (path if Path(path).is_file() else None)
    if not resolved:
        sys.exit(f"Could not find {name}: '{path}'. Provide its full path using the corresponding option.")
    return resolved


def collect_class_peptides(df):
    """Return class-I peptides, class-II peptides, and unsupported peptide lengths."""
    class1, class2, skipped = set(), set(), set()
    for col in PEPTIDE_COLUMNS.values():
        if col not in df.columns:
            sys.exit(f"Input file is missing the required column: '{col}'")
        for pep in df[col].dropna().astype(str).str.strip().str.upper():
            if not pep or not set(pep) <= VALID_AA:
                continue
            L = len(pep)
            if L in CLASS_I_LENGTHS:
                class1.add(pep)
            elif L in CLASS_II_LENGTHS:
                class2.add(pep)
            else:
                skipped.add(L)
    return sorted(class1), sorted(class2), skipped


# ---------------------------------------------------------------- multi-allele -xls parser

def parse_xls(path, tool, requested_alleles):
    """Parse NetMHCpan/NetMHCIIpan -xls output for one or more alleles.
    Columns are renamed to {tool}_{allele}_{name}. The header is identified by
    the 'Peptide' column; an optional allele-name row may appear immediately above it."""
    raw = [l.rstrip("\n") for l in Path(path).read_text().splitlines()]

    header = header_idx = None
    for i, line in enumerate(raw):
        toks = line.split("\t")
        if "Peptide" in toks:
            header, header_idx = toks, i
            break
    if header is None:
        first = "\n".join(raw[:5])
        sys.exit(f"Could not find a header containing the 'Peptide' column in {tool} output.\n"
                 f"File: {path}\nFirst lines:\n{first}")

    pep_idx = header.index("Peptide")

    # Is there an allele-name row immediately above the header?
    prev = raw[header_idx - 1].split("\t") if header_idx >= 1 else []
    allele_pos = [i for i, t in enumerate(prev)
                  if t.strip() and t.strip().lower() not in SKIP_COLS and t.strip().lower() not in LEAD_COLS]
    if allele_pos and "Peptide" not in prev:
        block_starts = allele_pos
        alleles = [prev[i].strip() for i in block_starts]
    else:
        # No allele-name row: assume a single output block and a single requested allele.
        if len(requested_alleles) != 1:
            first = "\n".join(raw[:5])
            sys.exit(f"The {tool} output has no allele-name row, but multiple alleles were requested: "
                     f"{requested_alleles}. First lines:\n{first}")
        start = next((j for j, h in enumerate(header) if h.strip().lower() not in LEAD_COLS), pep_idx + 1)
        block_starts = [start]
        alleles = [requested_alleles[0]]

    records = []
    for line in raw[header_idx + 1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) <= pep_idx:
            continue
        rec = {"Peptide": fields[pep_idx]}
        for k, start in enumerate(block_starts):
            end = block_starts[k + 1] if k + 1 < len(block_starts) else len(header)
            allele = alleles[k]
            for j in range(start, min(end, len(fields))):
                name = header[j].strip()
                if not name or name.lower() in SKIP_COLS or name == "Peptide":
                    continue
                rec[f"{tool}_{allele}_{name}"] = pd.to_numeric(fields[j], errors="coerce")
        records.append(rec)

    return pd.DataFrame(records).drop_duplicates("Peptide").set_index("Peptide")


# ---------------------------------------------------------------- predictor execution

def run_netmhcpan(exe, peptides, alleles, use_ba, workdir, tag, tool):
    pep_file = workdir / f"pep_{tool}_{tag}.txt"
    xls_file = workdir / f"{tool}_{tag}.xls"
    pep_file.write_text("\n".join(peptides) + "\n")
    cmd = [exe, "-p", str(pep_file), "-a", ",".join(alleles), "-xls", "-xlsfile", str(xls_file)]
    if use_ba:
        cmd.append("-BA")
    print("Running:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not xls_file.exists():
        sys.exit(f"NetMHCpan failed:\n{res.stderr or res.stdout}")
    return parse_xls(xls_file, tool, alleles)


def run_netmhc2pan(exe, peptides, alleles, use_ba, workdir, tag):
    pep_file = workdir / f"pep_net2_{tag}.txt"
    xls_file = workdir / f"net2_{tag}.xls"
    pep_file.write_text("\n".join(peptides) + "\n")
    cmd = [exe, "-f", str(pep_file), "-inptype", "1", "-a", ",".join(alleles),
           "-xls", "-xlsfile", str(xls_file)]
    if use_ba:
        cmd.append("-BA")
    print("Running:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not xls_file.exists():
        sys.exit(f"NetMHCIIpan failed:\n{res.stderr or res.stdout}")
    return parse_xls(xls_file, "netmhc2pan", alleles)


def run_mixmhcpred(exe, peptides, alleles, workdir, tag):
    """Run MixMHCpred. Output includes peptide, best-allele fields, scores and ranks;
    lines beginning with # are treated as comments."""
    in_file = workdir / f"pep_mix_{tag}.txt"
    out_file = workdir / f"mix_{tag}.txt"
    in_file.write_text("\n".join(peptides) + "\n")
    cmd = [exe, "-i", str(in_file), "-o", str(out_file), "-a", ",".join(alleles)]
    print("Running:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_file.exists():
        sys.exit(f"MixMHCpred failed:\n{res.stderr or res.stdout}")

    out = pd.read_csv(out_file, sep="\t", comment="#").drop_duplicates("Peptide").set_index("Peptide")
    rename = {}
    for c in out.columns:
        if c in ("Score_bestAllele", "%Rank_bestAllele", "BestAllele"):
            rename[c] = f"mixmhcpred_best_{c.replace('_bestAllele', '').replace('BestAllele', 'allele')}"
        elif c.startswith("Score_"):
            rename[c] = f"mixmhcpred_{c[len('Score_'):]}_Score"
        elif c.startswith("%Rank_"):
            rename[c] = f"mixmhcpred_{c[len('%Rank_'):]}_%Rank"
    out = out.rename(columns=rename)
    return out[[c for c in out.columns if c.startswith("mixmhcpred_")]]


def _strip_allele_suffix(name, allele):
    """Remove the _<allele> suffix from a column name, including unhyphenated forms."""
    for a in (allele, allele.replace("-", "")):
        suf = "_" + a
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def run_mixmhc2pred(exe, peptides, allele, with_context, workdir, tag):
    """Run MixMHC2pred. Prefix output columns with the requested allele explicitly,
    rather than relying on the allele naming convention used in its output."""
    in_file = workdir / f"pep_mix2_{tag}.txt"
    out_file = workdir / f"mix2_{tag}.txt"
    in_file.write_text("\n".join(peptides) + "\n")
    cmd = [exe, "-i", str(in_file), "-o", str(out_file), "-a", allele]
    if not with_context:
        cmd.append("--no_context")
    print("Running:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_file.exists():
        sys.exit(f"MixMHC2pred failed:\n{res.stderr or res.stdout}")

    out = pd.read_csv(out_file, sep=None, engine="python", comment="#")
    if "Peptide" not in out.columns:
        head = "\n".join(out_file.read_text().splitlines()[:5])
        sys.exit(f"MixMHC2pred output does not contain a 'Peptide' column. First lines:\n{head}")
    out["Peptide"] = out["Peptide"].astype(str).str.strip().str.upper()
    out = out.drop_duplicates("Peptide").set_index("Peptide")
    out = out.drop(columns=[c for c in out.columns if c.strip().lower() == "context"], errors="ignore")
    out.columns = [f"mixmhc2pred_{allele}_{_strip_allele_suffix(c, allele)}" for c in out.columns]
    return out


def run_hmmhc(exe, peptides, allele, workdir, tag):
    """Run hmmhc-predict. Expected CSV output includes peptide, -10logOdds,
    and percentile_rank (lower percentile rank indicates stronger predicted binding)."""
    in_file = workdir / f"pep_hmmhc_{tag}.txt"
    out_file = workdir / f"hmmhc_{tag}.csv"
    in_file.write_text("\n".join(peptides) + "\n")
    cmd = [exe, "--allele", allele, "--input", str(in_file), "--output", str(out_file)]
    print("Running:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_file.exists():
        sys.exit(f"hmmhc failed:\n{res.stderr or res.stdout}")
    out = pd.read_csv(out_file)
    pep_col = out.columns[0]
    out = out.drop_duplicates(pep_col).set_index(pep_col)
    out.columns = [f"hmmhc_{allele}_{c}" for c in out.columns]
    out.index.name = "Peptide"
    return out


# ---------------------------------------------------------------- prediction assembly

def build_predictions(class1, class2, args, workdir):
    frames = []
    use_ba = not args.no_ba

    if class1:
        exe = check_executable(args.netmhcpan, "NetMHCpan")
        frames.append(run_netmhcpan(exe, class1, NET_ALLELES, use_ba, workdir, "c1", "netmhcpan"))
        if not args.skip_mixmhcpred:
            exe = check_executable(args.mixmhcpred, "MixMHCpred")
            frames.append(run_mixmhcpred(exe, class1, MIX_ALLELES, workdir, "c1"))

    if class2:
        if not args.skip_netmhc2pan:
            exe = check_executable(args.netmhc2pan, "NetMHCIIpan")
            frames.append(run_netmhc2pan(exe, class2, NET2_ALLELES, use_ba, workdir, "c2"))
        if not args.skip_hmmhc:
            exe = check_executable(args.hmmhc, "hmmhc-predict")
            for al in HMMHC_ALLELES:
                frames.append(run_hmmhc(exe, class2, al, workdir, f"c2_{al}"))
        if args.mixmhc2pred and not args.skip_mixmhc2pred:
            exe = check_executable(args.mixmhc2pred, "MixMHC2pred")
            for al in MIX2_ALLELES:
                frames.append(run_mixmhc2pred(exe, class2, al, args.mixmhc2pred_context, workdir, f"c2_{al}"))

    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        return None
    return pd.concat(frames, axis=1)  # Concatenate columns, aligned by peptide.


def allele_subframe(preds, group):
    """Select columns for an allele group and remove the tool-specific allele name:
    {tool}_{tool_allele}_{metric} -> {tool}_{metric}."""
    cols = {}
    for tool, tool_allele in group:
        pref = f"{tool}_{tool_allele}_"
        for c in preds.columns:
            if c.startswith(pref):
                cols[c] = f"{tool}_{c[len(pref):]}"
    if not cols:
        return None
    return preds[list(cols)].rename(columns=cols)


def attach(df, preds):
    """Append prediction columns by matching WT/MT peptide sequences."""
    out = df
    for prefix, col in PEPTIDE_COLUMNS.items():
        key = df[col].astype(str).str.strip().str.upper()
        block = preds.reindex(key.values)
        block.columns = [f"{prefix}_{c}" for c in block.columns]
        block.index = df.index
        out = pd.concat([out, block], axis=1)
    return out


def process_file(in_path, args):
    df = pd.read_csv(in_path, sep="\t")
    has_index = not isinstance(df.index, pd.RangeIndex)
    print(f"\n=== {in_path}: {len(df)} rows ===", file=sys.stderr)

    # Remove prediction columns from previous runs; never modify the input file.
    old = [c for c in df.columns if OLD_COL_RE.match(c)]
    base = df.drop(columns=old) if old else df

    class1, class2, skipped = collect_class_peptides(base)
    if class1:
        print(f"  Class I (8–12 aa): {len(class1)} peptides -> NetMHCpan + MixMHCpred {NET_ALLELES}",
              file=sys.stderr)
    if class2:
        print(f"  Class II (15 aa): {len(class2)} peptides -> NetMHCIIpan + hmmhc {NET2_ALLELES}",
              file=sys.stderr)
    if skipped:
        print(f"  Skipped unsupported peptide lengths: {sorted(skipped)}", file=sys.stderr)

    if not class1 and not class2:
        print("  No peptides with supported lengths; no output files were created.", file=sys.stderr)
        return

    with tempfile.TemporaryDirectory() as tmp:
        preds = build_predictions(class1, class2, args, Path(tmp))
    if preds is None:
        print("  No predictions were generated; no output files were created.", file=sys.stderr)
        return

    # Track which predictors were enabled in this run.
    enabled = {"netmhcpan"}
    if not args.skip_mixmhcpred:
        enabled.add("mixmhcpred")
    if not args.skip_netmhc2pan:
        enabled.add("netmhc2pan")
    if not args.skip_hmmhc:
        enabled.add("hmmhc")
    if args.mixmhc2pred and not args.skip_mixmhc2pred:
        enabled.add("mixmhc2pred")

    # Write a separate output file for each allele.
    for allele, group in OUTPUT_GROUPS.items():
        sub = allele_subframe(preds, group)
        if sub is None or sub.empty:
            continue
        # Warn only when an enabled predictor did not produce expected columns.
        missing = [tool for tool, tool_allele in group
                   if tool in enabled
                   and not any(c.startswith(f"{tool}_{tool_allele}_") for c in preds.columns)]
        if missing:
            print(f"  WARNING [{allele}]: no columns returned by {missing}; "
                  f"check the allele naming used by this predictor.", file=sys.stderr)
        out_df = attach(base.copy(), sub)
        out_path = in_path.with_name(f"{in_path.stem}.{allele}{in_path.suffix}")
        out_df.to_csv(out_path, sep="\t", index=has_index)
        print(f"  -> {out_path.name} ({out_df.shape[1]} columns)", file=sys.stderr)


def main():
    args = parse_args()
    for raw in args.input:
        process_file(Path(raw).expanduser(), args)


if __name__ == "__main__":
    main()
