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
- [Examples](#examples)
- [Use case scenario data](#use-case-scenario-data)
- [Contributing](#contributing)  
- [License](#license)  

---

## Features

1. Quality Control
2. Robust analysis and meta-analysis of GWAS
3. Standard GWAS meta-analysis
4. Bayesian meta-analysis
5. Meta-analysis with summary statistics imputation
6. Functional enrichment analysis
 
### 1. Discrete Phenotypes

A tab-delimited file with the following columns:

| SNP  | CHR | BP        | aa1 | ab1 | bb1 | aa0 | ab0 | bb0 |
|------|-----|-----------|-----|-----|-----|-----|-----|-----|
| rs1  | 11  | 84095494  | 0   | 23  | 244 | 9   | 83  | 178 |
| rs2  | 17  | 19683106  | 8   | 76  | 183 | 31  | 121 | 118 |
| rs3  | 4   | 68129844  | 4   | 58  | 205 | 13  | 93  | 164 |
| rs4  | 10  | 44481115  | 2   | 37  | 228 | 8   | 76  | 186 |
| rs5  | 4   | 68116450  | 21  | 109 | 137 | 43  | 129 | 98  |

 

### 2. Continuous Phenotypes

A tab-delimited file with the following columns:

| SNP        | CHR | BP  | xaa      | sdaa   | naa  | xab       | sdab   | nab  | xbb      | sdbb   | nbb  |
|------------|-----|-----|----------|--------|------|-----------|--------|------|----------|--------|------|
| rs1000000  | 1   | 100 | -0.02762 | 0.9914 | 283  | -0.01343  | 0.9956 | 1866 | 0.0106   | 1.003  | 3102 |
| rs1000000  | 1   | 100 | -0.1341  | 1.053  | 45   | -0.01535  | 0.9676 | 262  | 0.01574  | 1.013  | 514  |
| rs1000000  | 1   | 100 | -0.3615  | 1.053  | 10   | 0.06759   | 0.9186 | 93   | -0.01571 | 1.040  | 170  |
| rs1000000  | 1   | 100 | -0.02187 | 0.9836 | 172  | -0.002833 | 1.001  | 1320 | 0.00291  | 1.001  | 2578 |
| rs10000010 | 1   | 200 | -0.0109  | 1.001  | 1354 | 0.006966  | 0.9842 | 2628 | -0.002799| 1.032  | 126  |

  
### 3. BETA/SE Format

A tab-delimited file with these columns:

| SNP        | CHR | BP        | A1 | A2 | BETA         | SE           |
|------------|-----|-----------|----|----|--------------|--------------|
| rs374450   | 1   | 12305285  | A  | G  | -0.117157114 | 0.170367578  |
| rs10489599 | 1   | 16585817  | A  | G  | 0.111312914  | 0.125714433  |
| rs9725656  | 1   | 18460525  | T  | C  | -0.038092709 | 0.123232082  |
| rs1106839  | 1   | 19420391  | A  | G  | 0.004437813  | 0.124570835  |
| rs2274001  | 1   | 19464772  | T  | C  | 0.012144373  | 0.124526074  |


  
- **Robust methods**: MIN, MAX, MERT, or FAST (MinP, Cauchy, CMC, MCM Combination tests).  
- **Bayesian meta-analysis**  
- **Summary statistics imputation** with LD statistics as reference panels (`pred_ld.py`).  
 
---

## Requirements

 
 - A Linux-based OS (e.g. Ubuntu 20.04 LTS)
- Python 3.10+  
- Requirements specified in requirements.txt
- `pyrama_beta_se_meta` binary or script (for beta/SE meta-analysis), compile :
```
g++ -std=c++17 -O3 -pthread -o pyrama_beta_se_meta PYRAMA_beta_SE_meta.cpp
```
- 'ref' folder containing LD information of each reference panel in parquet file format (available for download here : ... ). After downloading, extract to the working folder of PYRAMA. The Pheno Scanner LD reference panel is available upon request from the Pheno Scanner database (http://www.phenoscanner.medschl.cam.ac.uk) 


---

## Installation

1. Clone the repository:  
   ```bash
   git clone https://github.com/pbagos/PYRAMA.git
   cd PYRAMA
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
   Compile with the following command: 
   ```
   g++ -std=c++17 -O3 -pthread -o pyrama_beta_se_meta PYRAMA_beta_SE_meta.cpp
   ```

5. If you wish to download the LD reference panels to conduct summary statistics imputation with PRED-LD, visit this link: [link here...]

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
| `--robust_method`             | MIN, MAX, MERT, or FAST (MinP, Cauchy, MCM and CMC combination tests)                       |
| `--type_of_effect`            | FIXED or RANDOM                                                                             |
| `--approximate_max`           | YES or NO                                                                                   |
| `--bayesian_meta`             | YES to run Bayesian meta-analysis (default NO). Available only for Discrete Phenotype Input |
| `--imputation`                | Summary statistics imputation of missing SNPs  (only for BETA/SE input)                     |
| `--r2threshold`               | R² threshold for  imputation                                                                |
| `--population`                | Population code (e.g. EUR, AFR) for imputation                                              |
| `--maf`                       | Minor allele frequency cutoff for imputation                                                |
| `--ref`                       | Reference panel identifier (e.g. TOP_LD, Pheno_Scanner, Hap_Map, all_panels[default])       |
| `--missing_threshold`         | Fraction of studies that must include a SNP (default 0.5). When set to 0, it performs  imputation in all given variants|
| `-n`, `--nthreads`            | Number of parallel threads for file I/O (default: 1)                                        |

#### Case 2: Direct BETA/SE (No Imputation)

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]   --o beta_se_meta.txt    
```

#### Case 2: BETA/SE + Imputation

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]   --o imputed_meta.txt    --imputation   --r2threshold 0.8   --population EUR   --maf 0.01   --ref all_panels   --missing_threshold 0.0    
```

#### Case 1: Discrete Counts (Standard, Fast, Bayesian)

- **Standard**  
  ```bash
  python pyrama.py     --i study1.txt  study2.txt [study3.txt ...]    --o standard_meta.txt     --inheritance_model DOMINANT     --effect_size_type CATT     --robust_method MIN     --type_of_effect FIXED  
  ```
- **Fast Robust**  
  ```bash
  python pyrama.py     --i study1.txt  study2.txt [study3.txt ...]    --o fast_meta.txt     --inheritance_model ADDITIVE     --effect_size_type OR     --robust_method FAST    
  ```
- **Bayesian**  
  ```bash
  python pyrama.py     --i study1.txt  study2.txt [study3.txt ...]     --o bayes_meta.txt     --bayesian_meta YES     --inheritance_model ADDITIVE     --effect_size_type OR     
  ```

#### Case 3: Continuous Phenotype

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]    --o continuous_meta.txt   --inheritance_model ADDITIVE   --robust_method MAX   --type_of_effect RANDOM
```
## Use case scenario data 

The replication report is available at analysis_code.html  file of this repository with all data and results.

## Contributing

Contributions welcome! Please:

1. Fork this repo  
2. Create a feature branch  
3. Open a pull request with tests/examples  

---

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.
