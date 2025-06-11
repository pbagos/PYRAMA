 # PYRAMA

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A Python tool for robust analysis and meta-analysis of Genome Wide Association Studies.

---

## Table of Contents

- [Features](#features)  
- [Requirements](#requirements)  
- [Installation](#installation)  
- [Usage](#usage)  
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

- **Flexible input formats**: beta/SE, discrete phenotypes and continuous phenotypes.  
- **Robust methods**: MIN, MAX, MERT, or FAST (MinP, Cauchy,CMC,MCM Combination tests).  
- **Bayesian meta-analysis**  
- **Summary statistics Imputation** with LD statistics as reference panels (`pred_ld.py`).  
 
---

## Requirements

 
 - A Linux-based OS (e.g. Ubuntu 20.04 LTS)
- Python 3.10+  
- [pandas](https://pandas.pydata.org/)  
- `meta_analysis`, `cont_meta_analysis`, `fast_robust_analysis`, `bayesian`, `bivariate`, `bivariate_gwas` modules (included in this package).  
- `pred_ld.py` (for imputation)  
- `pyrama_beta_se_meta` binary or script (for beta/SE meta-analysis)  
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

 
  
### GWAS Meta-Analysis  
```bash
python pyrama.py   --i study1.txt study2.txt [study3.txt ...]   --o output_results.txt   [options]
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
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]   --o beta_se_meta.txt    
```

#### Case 2: BETA/SE + Imputation

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]   --o imputed_meta.txt    --imputation   --r2threshold 0.8   --population EUR   --maf 0.01   --ref 1000G   --missing_threshold 0.0    
```

#### Case 1: Discrete Counts (Standard, Fast, Bayesian)

- **Standard**  
  ```bash
  python pyrama.py     --i counts1.tsv counts2.txt     --o standard_meta.txt     --inheritance_model DOMINANT     --effect_size_type CATT     --robust_method MIN     --type_of_effect FIXED  
  ```
- **Fast Robust**  
  ```bash
  python pyrama.py     --i counts1.tsv counts2.txt     --o fast_meta.txt     --inheritance_model ADDITIVE     --effect_size_type OR     --robust_method FAST    
  ```
- **Bayesian**  
  ```bash
  python pyrama.py     --i counts1.tsv counts2.txt     --o bayes_meta.txt     --bayesian_meta YES     --inheritance_model ADDITIVE     --effect_size_type OR     
  ```

#### Case 3: Continuous Phenotype

```bash
python pyrama.py   --i cont1.tsv cont2.tsv cont3.txt   --o continuous_meta.txt   --inheritance_model ADDITIVE   --robust_method MAX   --type_of_effect RANDOM
```
 

## Contributing

Contributions welcome! Please:

1. Fork this repo  
2. Create a feature branch  
3. Open a pull request with tests/examples  

---

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.
