# PYRAMA
 
---

 
## Introduction

The **PYRAMA** (Python Robust Analysis and Meta-Analysis) toolkit is designed to streamline large-scale meta-analyses of GWAS summary statistics. This script handles multiple inheritance models, effect-size types, and phenotypic inputs (discrete, continuous, and bivariate), and it supports optional Bayesian approaches, robust methods (including FAST robust analysis), and genotype imputation for missing data. The script dynamically detects the input columns in your GWAS summary files and dispatches the appropriate analysis pipeline.

---

## Features

- **Automatic Input Detection**  
  - Discrete phenotype (case/control counts)  
  - Continuous phenotype (means, standard deviations, and sample sizes)  
  - Bivariate meta-analysis (count-based or beta/SE-based)

- **Multiple Inheritance Models**  
  - Additive  
  - Recessive  
  - Dominant  

- **Effect-Size Types**  
  - Odds Ratio (OR)  
  - Cochran-Armitage Trend Test (CATT)

- **Robust Methods**  
  - MIN, MAX, or MERT  
  - FAST robust analysis  
  - Bayesian meta-analysis (optional)

- **Imputation**  
  - Detects SNPs with missing data  
  - Interfaces with an external imputation module (`imputation.py`)  
  - Supports specifying R² threshold, population, MAF, reference panel, and optional SNP list

- **Bivariate Analysis**  
  - Counts-based bivariate meta-analysis (`bivariate.biv_meta_analysis`)  
  - Beta/SE-based bivariate meta-analysis (`bivariate_gwas.beta_SE_meta`)

- **Flexible Command‐Line Interface**  
  - Customizable input/output file paths  
  - Specify inheritance models, robust methods, effect types, etc.  
  - Optional imputation and Bayesian workflow  

- **Logging & Error Handling**  
  - Graceful detection of missing or malformed columns  
  - Helpful error messages for reading/writing files  

---

## Dependencies & Requirements

- **Python Version:**  
  - Tested with Python 3.8, 3.9, and 3.10 (≥3.8 recommended)

- **Python Packages:**  
  - `pandas` (≥1.0.0)  
  - `dask` (for `dask.dataframe` import; ≥2021.4.0)  
  - `argparse` (Standard library)  
  - Custom modules (included in the PYRAMA repo or installed separately):  
    - `meta_analysis`  
    - `beta_se_meta2`  
    - `imputation`  
    - `cont_meta_analysis`  
    - `fast_robust_analysis`  
    - `bayesian`  
    - `bivariate`  
    - `bivariate_gwas`  

- **External Dependencies (optional):**  
  - A compiled binary/executable named `pyrama_beta_SE_meta` for Beta/SE‐based meta-analysis.  
  - An imputation pipeline script (`imputation.py`) that accepts parameters:  
    ```
    python3 imputation.py <input_file> <output_dir> <r2threshold> <population> <maf> <ref> [imp_list]
    ```

---

## Installation

1. **Clone or Download the Repository**  
   ```bash
   git clone https://github.com/your_username/PYRAMA.git
   cd PYRAMA
