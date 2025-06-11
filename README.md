 # PYRAMA

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A Python tool for robust analysis and meta-analysis of genome-wide association studies (GWAS).

---

## Table of Contents

- [Features](#features)  
- [Requirements](#requirements)  
- [Installation](#installation)  
- [Usage](#usage)  
  - [Filtering Low-Frequency SNPs](#filtering-low-frequency-snps)  
  - [Merging Multiple TSV Files](#merging-multiple-tsv-files)  
  - [GWAS Meta-Analysis (Main CLI)](#gwas-meta-analysis-main-cli)  
    - [Case 2: Direct BETA/SE (No Imputation)](#case-2-direct-betase-no-imputation)  
    - [Case 2: BETA/SE + Imputation](#case-2-betase--imputation)  
    - [Case 1: Discrete Counts (Standard, Fast, Bayesian)](#case-1-discrete-counts-standard-fast-bayesian)  
    - [Case 3: Continuous Phenotype](#case-3-continuous-phenotype)  
    - [Bivariate Meta-Analysis](#bivariate-meta-analysis)  
- [API Reference](#api-reference)  
- [Examples](#examples)  
- [Contributing](#contributing)  
- [License](#license)  

---

## Features

- **Flexible input formats**: allele counts, β/SE, continuous traits.  
- **Robust methods**: MIN, MAX, MERT, or FAST.  
- **Bayesian meta-analysis** support.  
- **Missing data imputation** via LD‐based prediction (`pred_ld.py`).  
- **Bivariate analysis** (allele counts or β/SE).  
- **Parallel file reading** with thread pooling.  

---

## Requirements

- Python 3.7+  
- [pandas](https://pandas.pydata.org/)  
- `meta_analysis`, `cont_meta_analysis`, `fast_robust_analysis`, `bayesian`, `bivariate`, `bivariate_gwas` modules (included in this package).  
- `pred_ld.py` (for imputation)  
- `pyrama_beta_se_meta` binary or script (for β/SE meta-analysis)  

---

## Installation

1. Clone the repository:  
   ```bash
   git clone https://github.com/yourusername/pyrama.git
   cd pyrama
   ```
2. (Optional) Create and activate a virtual environment:  
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
3. Install dependencies:  
   ```bash
   pip install -r requirements.txt
   ```
4. Ensure `pred_ld.py` is executable and on your `PATH`, and `pyrama_beta_se_meta` is compiled.

---

## Usage

All examples assume your script is called `pyrama.py`. Adjust accordingly.

### Filtering Low-Frequency SNPs

Use `process_snps()` to write all SNPs that occur less than the maximum frequency in a TSV:

```bash
python - <<'EOS'
from pyrama import process_snps
process_snps(
    input_file="all_snps.tsv",
    output_file="low_freq_snps.txt"
)
EOS
```

### Merging Multiple TSV Files

Programmatic merge of many studies into a single DataFrame, sorted by `SNP`:

```python
from pyrama import merge_input_files

df = merge_input_files(
    file_list=[
        "study1.tsv",
        "study2.tsv",
        "study3.tsv"
    ],
    max_workers=4
)
print(df.head())
```

### GWAS Meta-Analysis (Main CLI)

```bash
python pyrama.py   --i study1.tsv study2.tsv [study3.tsv ...]   --o output_results.tsv   [options]
```

Common options:

| Flag                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--inheritance_model`         | ADDITIVE, RECESSIVE, or DOMINANT                                                            |
| `--effect_size_type`          | OR (odds ratio) or CATT (Cochran–Armitage trend test)                                       |
| `--robust_method`             | MIN, MAX, MERT, or FAST                                                                     |
| `--type_of_effect`            | FIXED or RANDOM                                                                             |
| `--approximate_max`           | YES or NO                                                                                   |
| `--biv_ma`                    | YES to force bivariate meta-analysis (default NO)                                           |
| `--bayesian_meta`             | YES to run Bayesian meta-analysis (default NO)                                              |
| `--imputation`                | Enable imputation of missing SNPs                                                           |
| `--r2threshold`               | R² threshold for LD‐based imputation                                                        |
| `--population`                | Population code (e.g. EUR, AFR) for imputation                                              |
| `--maf`                       | Minimum allele frequency cutoff for imputation                                              |
| `--ref`                       | Reference panel identifier (e.g. 1000G)                                                     |
| `--missing_threshold`         | Fraction of studies that must include a SNP (default 0.5)                                   |
| `-n`, `--nthreads`            | Number of parallel threads for file I/O (default: 1)                                        |

#### Case 2: Direct BETA/SE (No Imputation)

```bash
python pyrama.py   --i studyA.tsv studyB.tsv   --o beta_se_meta.tsv   --inheritance_model ADDITIVE   --effect_size_type OR   --robust_method MERT   --type_of_effect FIXED   --approximate_max NO
```

#### Case 2: BETA/SE + Imputation

```bash
python pyrama.py   --i studyA.tsv studyB.tsv studyC.tsv   --o imputed_meta.tsv   --inheritance_model RECESSIVE   --effect_size_type OR   --robust_method FAST   --type_of_effect RANDOM   --approximate_max YES   --imputation   --r2threshold 0.8   --population EUR   --maf 0.01   --ref 1000G   --missing_threshold 0.0   --nthreads 4
```

#### Case 1: Discrete Counts (Standard, Fast, Bayesian)

- **Standard**  
  ```bash
  python pyrama.py     --i counts1.tsv counts2.tsv     --o standard_meta.tsv     --inheritance_model DOMINANT     --effect_size_type CATT     --robust_method MIN     --type_of_effect FIXED     --approximate_max NO
  ```
- **Fast Robust**  
  ```bash
  python pyrama.py     --i counts1.tsv counts2.tsv     --o fast_meta.tsv     --inheritance_model ADDITIVE     --effect_size_type OR     --robust_method FAST     --approximate_max NO
  ```
- **Bayesian**  
  ```bash
  python pyrama.py     --i counts1.tsv counts2.tsv     --o bayes_meta.tsv     --bayesian_meta YES     --inheritance_model ADDITIVE     --effect_size_type OR     --approximate_max YES
  ```

#### Case 3: Continuous Phenotype

```bash
python pyrama.py   --i cont1.tsv cont2.tsv cont3.tsv   --o continuous_meta.tsv   --inheritance_model ADDITIVE   --robust_method MAX   --type_of_effect RANDOM
```

### Bivariate Meta-Analysis

- **Allele Counts**  
  ```bash
  python pyrama.py     --i biv_counts1.tsv biv_counts2.tsv     --o biv_counts_results.tsv     --inheritance_model RECESSIVE     --effect_size_type OR     --robust_method MERT     --type_of_effect FIXED     --approximate_max NO     --biv_ma YES
  ```
- **β/SE**  
  ```bash
  python pyrama.py     --i biv_beta1.tsv biv_beta2.tsv     --o biv_beta_se_results.tsv     --biv_ma YES
  ```

---

## API Reference

- **`process_snps(input_file: str, output_file: str) -> None`**  
  Filter out SNPs occurring at the maximum frequency; writes lower‐frequency SNPs to `output_file`.

- **`merge_input_files(file_list: List[str], max_workers: int = 1) -> pd.DataFrame`**  
  Reads and concatenates multiple TSVs in parallel, sorts by `SNP`, and returns a DataFrame.

- **`gwas_meta_analysis(...) -> None`**  
  Main entry point for GWAS meta-analysis. See [Usage](#usage) for parameter details.

Refer to the docstrings in each module (`meta_analysis`, `cont_meta_analysis`, etc.) for advanced options.

---

## Examples

1. **Filter SNPs**  
   ```bash
   python - <<'EOS'
   from pyrama import process_snps
   process_snps("study_all_snps.tsv", "rare_snps.txt")
   EOS
   ```
2. **Merge Files & Inspect**  
   ```python
   from pyrama import merge_input_files
   df = merge_input_files(["s1.tsv","s2.tsv"], max_workers=2)
   print(df.shape)
   ```
3. **Run Full Meta-Analysis**  
   ```bash
   python pyrama.py      --i s1.tsv s2.tsv s3.tsv      --o final_meta.tsv      --inheritance_model ADDITIVE      --effect_size_type OR      --robust_method MAX      --type_of_effect RANDOM      --approximate_max YES      --imputation      --r2threshold 0.7      --population AFR      --maf 0.05      --ref 1000G      --missing_threshold 0.2      --nthreads 8
   ```

---

## Contributing

Contributions welcome! Please:

1. Fork this repo  
2. Create a feature branch  
3. Open a pull request with tests/examples  

---

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.
