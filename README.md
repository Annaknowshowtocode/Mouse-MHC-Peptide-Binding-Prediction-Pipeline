# Mouse MHC Peptide-Binding Prediction Pipeline

A command-line Python pipeline for running and consolidating **in silico peptide–MHC binding predictions** for mouse MHC alleles. It was developed for the research project **“Immunometabolic control of tumor antigen presentation”** (Tomas Paulenda, Artyomov).

The pipeline takes TSV tables containing wild-type (WT) and mutant (MT) peptide sequences, routes peptides to class-I or class-II prediction tools according to peptide length, runs the configured predictors, and writes allele-specific result tables while preserving the original input columns.

> **Scope and reproducibility note:** this repository contains the Python orchestration and result-processing code. NetMHCpan, MixMHCpred, NetMHCIIpan, hmmhc, and MixMHC2pred are external tools and are **not bundled**. Their installation, model data, platform support, and licensing/access conditions are managed separately. Therefore, cloning this repository alone is not sufficient to run predictions.

## What the pipeline does

- Reads one or more tab-separated (`.tsv`) input files.
- Extracts peptide sequences from `wt_short_peptide` and `mt_short_peptide`.
- Validates sequences against the standard amino-acid alphabet and removes duplicate peptides before prediction.
- Assigns peptides by length:
  - **MHC class I:** 8–12 amino acids; predicts against both **H-2-Kb** and **H-2-Db**.
  - **MHC class II:** 15 amino acids; predicts against **H-2-IAb**.
  - Other lengths are skipped and reported.
- Integrates outputs from multiple prediction tools:
  | MHC class | Tool | Role in this pipeline |
  |---|---|---|
  | I | NetMHCpan 4.2 | Binding prediction for H-2-Kb and H-2-Db |
  | I | MixMHCpred 3.0 | Complementary prediction for H-2-Kb and H-2-Db |
  | II | NetMHCIIpan 4.3 | Binding prediction for H-2-IAb |
  | II | hmmhc | Additional class-II prediction for H2-IAb |
  | II | MixMHC2pred | Optional additional class-II prediction for H-2-IAb |
- Handles tool-specific allele naming conventions when consolidating output columns.
- Matches predictions back to WT and MT peptide columns and retains the original input data.
- Writes a separate output TSV for each allele; input files are not overwritten.

## Repository contents

```text
mhc_prediction_pipeline/
├── predict_mhc.py       # Main pipeline
├── requirements.txt     # Python dependency
├── .gitignore           # Generated files and local environments
└── README.md            # This documentation
```

## Requirements

### Python

- Python 3.9 or newer
- `pandas`

Install the Python dependency in a virtual environment:

```bash
# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows, use WSL2 or a Linux environment supported by the external prediction tools. The external tools may have their own operating-system, compiler, model-data, and license requirements.

### External prediction tools

Install and configure the tools independently, following their official documentation and the access terms supplied by their developers:

- NetMHCpan 4.2
- MixMHCpred 3.0
- NetMHCIIpan 4.3
- `hmmhc-predict`
- MixMHC2pred (optional)

The executable files must be available on `PATH`, or their full paths must be passed as command-line arguments. Verify that each tool works independently before running this pipeline.

## Input format

Each input file must be a tab-separated table with at least these two columns:

| Column | Description |
|---|---|
| `wt_short_peptide` | Wild-type peptide sequence |
| `mt_short_peptide` | Mutant peptide sequence |

Other columns are retained in the output. Peptide sequences are normalized to uppercase and stripped of surrounding whitespace for prediction. Rows with missing or non-standard amino-acid sequences are not submitted to the predictors.

Example:

```tsv
sample_id	wt_short_peptide	mt_short_peptide
sample_001	SIINFEKL	SIINFEKM
sample_002	ACDEFGHIKLMNPQR	ACDEFGHIKLMNPQS
```

The sequences above illustrate the expected column structure only; use biologically valid peptide sequences for actual analyses.

## Running the pipeline

Run the script from the repository directory or provide its full path. Supply one or more input files with `-i`:

```bash
python predict_mhc.py \
  -i data/sample_1.tsv data/sample_2.tsv \
  --netmhcpan /path/to/netMHCpan \
  --mixmhcpred /path/to/MixMHCpred \
  --netmhc2pan /path/to/netMHCIIpan \
  --hmmhc /path/to/hmmhc-predict
```

If the executables are already on `PATH`, their command names can be used instead:

```bash
python predict_mhc.py \
  -i data/sample_1.tsv \
  --netmhcpan netMHCpan \
  --mixmhcpred MixMHCpred \
  --netmhc2pan netMHCIIpan \
  --hmmhc hmmhc-predict
```

Add MixMHC2pred when it is installed and configured:

```bash
python predict_mhc.py \
  -i data/sample_1.tsv \
  --mixmhc2pred /path/to/MixMHC2pred \
  --netmhcpan /path/to/netMHCpan \
  --mixmhcpred /path/to/MixMHCpred \
  --netmhc2pan /path/to/netMHCIIpan \
  --hmmhc /path/to/hmmhc-predict
```

The script requires `-i`; there is no machine-specific default input path. Use `python predict_mhc.py --help` to see all options.

### Optional flags

| Flag | Effect |
|---|---|
| `--no-ba` | Do not request binding-affinity output from NetMHCpan / NetMHCIIpan |
| `--skip-mixmhcpred` | Skip MixMHCpred (class I) |
| `--skip-netmhc2pan` | Skip NetMHCIIpan (class II) |
| `--skip-hmmhc` | Skip hmmhc (class II) |
| `--skip-mixmhc2pred` | Skip MixMHC2pred (class II) |
| `--mixmhc2pred-context` | Pass context to MixMHC2pred; default behavior is `--no_context` |

If a class-II predictor is skipped, the corresponding output columns will not be present. MixMHC2pred is optional and is only run when `--mixmhc2pred` is supplied.

## Output files

For each input file, the pipeline creates one output file per allele for which predictions were generated. The output is written **next to the input file**, and the original input remains unchanged.

For example, an input named `candidate_peptides.tsv` may produce:

```text
candidate_peptides.H-2-Kb.tsv
candidate_peptides.H-2-Db.tsv
candidate_peptides.H-2-IAb.tsv
```

Each output contains:

1. All original input columns.
2. WT prediction columns prefixed with `wt_`.
3. MT prediction columns prefixed with `mt_`.
4. Predictor-specific metrics, with allele-specific predictions placed in the corresponding allele file.

The exact metrics depend on the output provided by each external tool and the options used. Missing predictions remain empty after peptide-based matching. The pipeline does not interpret a score as a biological validation or experimental confirmation.

## How to read the results

1. **Start with the allele in the filename.** `H-2-Kb`, `H-2-Db`, and `H-2-IAb` are separate result sets; compare predictions within the relevant allele rather than mixing allele-specific columns.
2. **Find the original peptide columns.** `wt_short_peptide` and `mt_short_peptide` identify the sequences being evaluated.
3. **Compare WT and MT values.** Columns beginning with `wt_` refer to the wild-type peptide; columns beginning with `mt_` refer to the mutant peptide.
4. **Identify the predictor and metric from the column name.** For example, `wt_netmhcpan_*` is a NetMHCpan metric for the WT peptide, while `mt_mixmhcpred_*` is a MixMHCpred metric for the MT peptide. Exact suffixes are inherited from the external tool's output.
5. **Interpret each metric using that predictor's documentation.** Scores and percentile ranks are not interchangeable across tools. In particular, do not compare raw scores from different predictors as though they shared one scale.
6. **Treat predictions as prioritization evidence.** Binding predictions can help prioritize candidate peptides, but do not by themselves establish antigen processing, cell-surface presentation, T-cell recognition, or immunogenicity.

## Reproducibility and portability

To run this project on another computer:

1. Clone or download this repository.
2. Install Python and the Python dependencies.
3. Install the external prediction tools and their required model/data files according to their own instructions.
4. Make the executables available on `PATH`, or provide their full paths in the command.
5. Prepare TSV files with the required peptide columns.
6. Run the command with explicit input and executable paths.

The Python code no longer contains a personal absolute input-file path. However, **the full workflow is not self-contained or universally portable**: execution depends on separately installed external predictors, their compatible model data, operating-system support, and access permissions. Do not commit licensed binaries, private model files, or sensitive input data to this repository.

For a reproducible analysis record, save the tool versions, command used, input-file checksums, and relevant predictor settings alongside the results.

## Limitations

- Peptide class is assigned by sequence length using the rules implemented here (8–12 aa for class I and 15 aa for class II).
- The pipeline orchestrates prediction tools; it does not perform mass-spectrometry validation, immunogenicity prediction, or experimental confirmation.
- Tool-specific metrics retain their original meaning and scale; cross-tool score comparisons require care.
- Results depend on the installed predictor versions, model data, allele support, and tool configuration.
- The pipeline expects the input column names documented above.

## Citation and attribution

If you use this pipeline in academic work, cite the original prediction tools according to their developers' recommendations and acknowledge the research project context:

**Project:** *Immunometabolic control of tumor antigen presentation*  
**Researchers provided for project attribution:** Tomas Paulenda; Artyomov.

Add the repository URL and the exact tool versions used when describing an analysis based on this pipeline.
