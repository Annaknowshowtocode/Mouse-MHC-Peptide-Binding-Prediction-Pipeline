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
├── predict_mhc.py       
├── requirements.txt    
├── .gitignore           
└── README.md           
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

### External prediction tools

Install and configure the tools independently, following their official documentation and the access terms supplied by their developers:

- NetMHCpan 4.2
- MixMHCpred 3.0
- NetMHCIIpan 4.3
- `hmmhc-predict`
- MixMHC2pred 

## Input format

Each input file must be a tab-separated table with at least these two columns:

| Column | Description |
|---|---|
| `wt_short_peptide` | Wild-type peptide sequence |
| `mt_short_peptide` | Mutant peptide sequence |

Example:

```tsv
sample_id	wt_short_peptide	mt_short_peptide
sample_001	SIINFEKL	SIINFEKM
sample_002	ACDEFGHIKLMNPQR	ACDEFGHIKLMNPQS
```

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

Use `python predict_mhc.py --help` to see all options.

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
