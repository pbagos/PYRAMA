 # PYRAMA

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A Python tool for robust analysis and meta-analysis of Genome Wide Association Studies.

---
 
## 📑 Table of Contents

- [Features](#features)
- [Input File Formats](#input-file-formats)
  - [Discrete Phenotypes](#discrete-phenotypes)
  - [Continuous Phenotypes](#continuous-phenotypes)
  - [BETA/SE Format](#betase-format)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Quality Control](#quality-control)
  - [GWAS Meta-Analysis](#gwas-meta-analysis)
    - [BETA/SE Format](#case-direct-betase-no-imputation)
      - [Without Imputation](#case-direct-betase-no-imputation)
      - [With Imputation](#case-betase--imputation)
    - [Discrete Phenotypes](#case-discrete-phenotypes-standard-fast-bayesian)
      - [Standard Meta-Analysis](#case-discrete-phenotypes-standard-fast-bayesian)
      - [Fast Robust Analysis](#case-discrete-phenotypes-standard-fast-bayesian)
      - [Bayesian Meta-Analysis](#case-discrete-phenotypes-standard-fast-bayesian)
    - [Continuous Phenotypes](#case-continuous-phenotypes)
- [Use Case Scenario Data](#use-case-scenario-data)
- [Contributing](#contributing)
- [License](#license)

## Features

1. Quality Control 
2. Robust analysis and meta-analysis of GWAS (For Discrete and Continuous phenotypes inputs)
3. Standard GWAS meta-analysis 
4. Bayesian meta-analysis (For Discrete and Continuous phenotypes inputs)
5. Meta-analysis with summary statistics imputation (only for BETA/SE input)
6. Functional enrichment analysis 

### File inputs examples 
#### 1. Discrete Phenotypes

A tab-delimited file with the following columns:

| SNP  | CHR | BP        | aa1 | ab1 | bb1 | aa0 | ab0 | bb0 |
|------|-----|-----------|-----|-----|-----|-----|-----|-----|
| rs1  | 11  | 84095494  | 0   | 23  | 244 | 9   | 83  | 178 |
| rs2  | 17  | 19683106  | 8   | 76  | 183 | 31  | 121 | 118 |
| rs3  | 4   | 68129844  | 4   | 58  | 205 | 13  | 93  | 164 |
| rs4  | 10  | 44481115  | 2   | 37  | 228 | 8   | 76  | 186 |
| rs5  | 4   | 68116450  | 21  | 109 | 137 | 43  | 129 | 98  |

 | Column | Description                                                              |
| ------ | ------------------------------------------------------------------------ |
| `SNP`  | SNP identifier (e.g., rsID).                                             |
| `CHR`  | Chromosome number where the SNP is located.                              |
| `BP`   | Base-pair position of the SNP on the chromosome (genomic coordinate).    |
| `aa1`  | Count of individuals with the homozygous reference genotype in cases.    |
| `ab1`  | Count of individuals with the heterozygous genotype in cases.            |
| `bb1`  | Count of individuals with the homozygous alternate genotype in cases.    |
| `aa0`  | Count of individuals with the homozygous reference genotype in controls. |
| `ab0`  | Count of individuals with the heterozygous genotype in controls.         |
| `bb0`  | Count of individuals with the homozygous alternate genotype in controls. |


#### 2. Continuous Phenotypes

A tab-delimited file with the following columns:

| SNP        | CHR | BP  | xaa      | sdaa   | naa  | xab       | sdab   | nab  | xbb      | sdbb   | nbb  |
|------------|-----|-----|----------|--------|------|-----------|--------|------|----------|--------|------|
| rs1000000  | 1   | 100 | -0.02762 | 0.9914 | 283  | -0.01343  | 0.9956 | 1866 | 0.0106   | 1.003  | 3102 |
| rs1000000  | 1   | 100 | -0.1341  | 1.053  | 45   | -0.01535  | 0.9676 | 262  | 0.01574  | 1.013  | 514  |
| rs1000000  | 1   | 100 | -0.3615  | 1.053  | 10   | 0.06759   | 0.9186 | 93   | -0.01571 | 1.040  | 170  |
| rs1000000  | 1   | 100 | -0.02187 | 0.9836 | 172  | -0.002833 | 1.001  | 1320 | 0.00291  | 1.001  | 2578 |
| rs10000010 | 1   | 200 | -0.0109  | 1.001  | 1354 | 0.006966  | 0.9842 | 2628 | -0.002799| 1.032  | 126  |

| Column | Description                                                                  |
| ------ | ---------------------------------------------------------------------------- |
| `SNP`  | SNP identifier (e.g., rsID).                                                 |
| `CHR`  | Chromosome number where the SNP is located.                                  |
| `BP`   | Base-pair position of the SNP on the chromosome.                             |
| `xaa`  | Mean phenotype value for individuals with the homozygous reference genotype. |
| `sdaa` | Standard deviation of the phenotype in the homozygous reference group.       |
| `naa`  | Number of individuals with the homozygous reference genotype.                |
| `xab`  | Mean phenotype value for individuals with the heterozygous genotype.         |
| `sdab` | Standard deviation of the phenotype in the heterozygous group.               |
| `nab`  | Number of individuals with the heterozygous genotype.                        |
| `xbb`  | Mean phenotype value for individuals with the homozygous alternate genotype. |
| `sdbb` | Standard deviation of the phenotype in the homozygous alternate group.       |
| `nbb`  | Number of individuals with the homozygous alternate genotype.                |
  
#### 3. BETA/SE Format

A tab-delimited file with these columns:

| SNP        | CHR | BP        | A1 | A2 | BETA         | SE           |
|------------|-----|-----------|----|----|--------------|--------------|
| rs374450   | 1   | 12305285  | A  | G  | -0.117157114 | 0.170367578  |
| rs10489599 | 1   | 16585817  | A  | G  | 0.111312914  | 0.125714433  |
| rs9725656  | 1   | 18460525  | T  | C  | -0.038092709 | 0.123232082  |
| rs1106839  | 1   | 19420391  | A  | G  | 0.004437813  | 0.124570835  |
| rs2274001  | 1   | 19464772  | T  | C  | 0.012144373  | 0.124526074  |

| Column | Description                                                          |
| ------ | -------------------------------------------------------------------- |
| `SNP`  | SNP identifier (e.g., rsID)                                         |
| `CHR`  | Chromosome number where the SNP is located                          |
| `BP`   | Base-pair position of the SNP on the chromosome                     |
| `A1`   | Effect allele (the allele associated with the reported effect size) |
| `A2`   | Non-effect allele (the alternative allele)                          |
| `BETA` | Estimated effect size (OR/log(OR)) |
| `SE`   | Standard error of the estimated effect size                         |
 
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

### Quality control
Before conducting an analysis or a meta-analysis, we recommend that users execute the quality control script that accompanies PYRAMA. This script performs essential preprocessing steps, including the removal of problematic rows from the input studies, verification of allele order (harmonization) and allele consistency across all datasets, and generation of a comprehensive quality control report. Additionally, it calculates the number of shared variants for each combination of studies, providing users with a clear overview of the variant over-lap among the imported studies. 

```bash
python quality_control .py  --input_files study1.txt study2.txt [study3.txt ...]    --output final_merged_gwas.txt  [--skip_harm (optional)]
```
 
| Flag                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--input_files` |	List of input GWAS files (tab-delimited). Each file should contain columns: SNP, CHR, BP, A1, A2, BETA, and SE. You must provide at least one input file.|
| `--output` | Output filename for the final merged and harmonized GWAS data. Will be written in tab-separated format.|
| `--skip_harm` | If set, skips the allele harmonization and filtering steps involving A1/A2. Use this when harmonization is not needed or when input lacks allele information.|
 
### GWAS Meta-Analysis  
```bash
python pyrama.py  --i study1.txt study2.txt [study3.txt ...]  --o output_results.txt   [options]
```

Common options:

| Flag                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--i`, `--input`, `--input_file`|  Paths to input data files (space-separated when multiple files are inserted)             |
| `--o`, `--output`, `--output_file`|  Path to output files                                                                    |
| `--inheritance_model`         | ADDITIVE, RECESSIVE, or DOMINANT                                                            |
| `--effect_size_type`          | OR (odds ratio) or CATT (Cochran–Armitage trend test)                                       |
| `--robust_method`             | MIN, MAX, MERT, or FAST (MinP, Cauchy, MCM and CMC combination tests)                       |
| `--type_of_effect`            | FIXED or RANDOM                                                                             |
| `--bayesian_meta`             | YES to run Bayesian meta-analysis (default NO). Available only for Discrete Phenotype Input |
| `--imputation`                | Summary statistics imputation of missing SNPs  (only for BETA/SE input)                     |
| `--r2threshold`               | R² threshold for  imputation                                                                |
| `--population`                | Population code (EUR, AFR,EAS and SAS) for imputation                                              |
| `--maf`                       | Minor allele frequency cutoff for imputation                                                |
| `--ref`                       | Reference panel identifier (e.g. TOP_LD, Pheno_Scanner, Hap_Map, all_panels[default])       |
| `--missing_threshold`         | Fraction of studies that must include a SNP (default 0.5). When set to 0, it performs  imputation in all given variants|
| `-n`, `--nthreads`            | Number of parallel threads for file I/O (default: 1)                                        |
| `--het_est`                   | Heterogeneity estimator: 'DL' (DerSimonian-Laird) [default], "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman) |

#### Case: Direct BETA/SE (No Imputation)

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]   --o beta_se_meta.txt    
```

#### Case: BETA/SE + Imputation

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]   --o imputed_meta.txt    --imputation   --r2threshold 0.8   --population EUR   --maf 0.01   --ref all_panels   --missing_threshold 0.0    
```
##### Output columns

| Column     | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
| `SNP`  | SNP identifier (e.g., rsID).                                                 |
| `CHR`  | Chromosome number where the SNP is located                                  |
| `BP`   | Base-pair position of the SNP on the chromosome                             |
| `N`          | Number of valid studies for this variant|
| `P`          | P-value of Random Effects meta-analysis|
| `SE`         | Standard error|
| `BETA`       | Random Effects estimate |
| `I2`         | Heterogeneity index (I² statistic) in Random Effects meta-analysis|
| `pQ`         | P-value for heterogeneity (often from Cochran's Q test)|
| `BETA(FE)`   | Fixed Effect estimate|
| `P(FE)`      | P-value for the Fixed Effect meta-analysis|




#### Case: Discrete Phenotypes (Standard, Fast, Bayesian)

| Flag                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--i`, `--input`, `--input_file`|  Paths to input data files (space-separated when multiple files are inserted)             |
| `--o`, `--output`, `--output_file`|  Path to output files                                                                    |
| `--inheritance_model`         | ADDITIVE, RECESSIVE, DOMINANT or ALL  (to enable robust methods)                             |
| `--effect_size_type`          | OR (Odds Ratio) or CATT (Cochran–Armitage trend test)                                       |
| `--robust_method`             | MIN (MIN2), MAX, MERT, or FAST (MinP, Cauchy, MCM and CMC combination tests)                       |
| `--type_of_effect`            | FIXED or RANDOM                                                                             |
| `--bayesian_meta`             | YES to run Bayesian meta-analysis (default NO). Available  for Discrete and Continuous Phenotypes Inputs |
| `-n`, `--nthreads`            | Number of parallel threads for file I/O (default: 1)                                        |



- **Standard**  
  ```bash
  python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]    --o standard_meta.txt   --inheritance_model ADDITIVE   --effect_size_type OR   --type_of_effect FIXED  
  ```

- **MIN2/MAX/MERT robust methods**  
  ```bash
  python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]    --o standard_meta.txt   --inheritance_model ALL   --effect_size_type OR     --robust_method MIN --type_of_effect RANDOM  
  ```


##### Output columns

| Column     | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
| `SNP`  | SNP identifier (e.g., rsID).                                                 |
| `CHR`  | Chromosome number where the SNP is located                                  |
| `BP`   | Base-pair position of the SNP on the chromosome                             |
| `P`    | P-value of the analysis/meta-analysis|
| `CI_lower`| Lower confidence interval|
| `CI_upper`| Upper confidence interval|
| `Weighted_Mean`|Analysis estimate|



  
- **Fast Robust**  
  ```bash
  python pyrama.py     --i study1.txt  study2.txt [study3.txt ...]    --o fast_meta.txt     --inheritance_model ALL     --effect_size_type OR     --robust_method FAST   --het_est DL
  ```


| Flag                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--i`, `--input`, `--input_file`|  Paths to input data files (space-separated when multiple files are inserted)             |
| `--o`, `--output`, `--output_file`|  Path to output files                                                                    |
| `--inheritance_model`         | ADDITIVE, RECESSIVE, or DOMINANT                                                            |
| `--effect_size_type`          | OR (Odds Ratio) or CATT (Cochran–Armitage trend test)                                       |
| `--robust_method`             | MIN, MAX, MERT, or FAST (MinP, Cauchy, MCM and CMC combination tests)                       |
| `--type_of_effect`            | FIXED or RANDOM                                                                             |
| `--bayesian_meta`             | YES to run Bayesian meta-analysis (default NO). Available  for Discrete and Continuous Phenotypes Inputs |
| `-n`, `--nthreads`            | Number of parallel threads for file I/O (default: 1)                                        |
| `--het_est`                   | Heterogeneity estimator: 'DL' (DerSimonian-Laird) [default], "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman) |

##### Output columns

| Column     | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
| `SNP`  | SNP identifier (e.g., rsID).                                                 |
| `CHR`  | Chromosome number where the SNP is located                                  |
| `BP`   | Base-pair position of the SNP on the chromosome                             |
| `Z_Dom`    | Z-score value of the dominant model|
| `Z_Add`    | Z-score value of the additive (allelic) model|
| `Z_Rec`    | Z-score value of the recessive model|
| `P_Dom`    | p-value from the Dominant model|
| `P_Add`    | p-value from the Additive (Allelic) model|
| `P_Rec`    | p-value from the Recessive model|
| `P_MinP`    | p-value from the MinP combination test|
| `P_CCT`    | p-value from the  Cauchy combination test (CCT)|
| `P_CMC`    | p-value from the  CMC combination test|
| `P_MCM`    | p-value from the  MCM combination test|
| `P_Stouffer` | p-value from the Stouffer's method for dependent tests combination test|




  
- **Bayesian**  
  ```bash
  python pyrama.py     --i study1.txt  study2.txt [study3.txt ...]     --o bayes_meta.txt     --bayesian_meta YES     --inheritance_model ADDITIVE     --effect_size_type OR     
  ```
##### Output columns

| Column     | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
| `SNP`  | SNP identifier (e.g., rsID).                                                 |
| `CHR`  | Chromosome number where the SNP is located                                  |
| `BP`   | Base-pair position of the SNP on the chromosome                             |
| `P`    | P-value of the analysis/meta-analysis|
| `Z`    | Z-score value of the analysis/meta-analysis|
| `CI_low`| Lower confidence interval|
| `CI_upp`| Upper confidence interval|
| `E_m`    | Posterior expectation for the population parameter μ |
| `E_tau_square`    | Posterior expectation of the between study variability τ^2|
| `V_mu`    | Variance of the posterior distribution for the population parameter μ |
| `V_tau_square`    | Posterior variance of the between study variability τ^2|
 



#### Case: Continuous Phenotypes

```bash
python pyrama.py   --i study1.txt  study2.txt [study3.txt ...]    --o continuous_meta.txt   --inheritance_model ADDITIVE   --robust_method MAX   --type_of_effect RANDOM
```
| Flag                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--i`, `--input`, `--input_file`|  Paths to input data files (space-separated when multiple files are inserted)             |
| `--o`, `--output`, `--output_file`|  Path to output files                                                                    |
| `--inheritance_model`         | ADDITIVE, RECESSIVE, or DOMINANT                                                            |
| `--effect_size_type`          | OR (Odds Ratio) or CATT (Cochran–Armitage trend test)                                       |
| `--robust_method`             | MIN, MAX, MERT, or FAST (MinP, Cauchy, MCM and CMC combination tests)                       |
| `--type_of_effect`            | FIXED or RANDOM                                                                             |
| `--bayesian_meta`             | YES to run Bayesian meta-analysis (default NO). Available  for Discrete and Continuous Phenotypes Inputs |
| `-n`, `--nthreads`            | Number of parallel threads for file I/O (default: 1)                                        |
| `--het_est`                   | Heterogeneity estimator: 'DL' (DerSimonian-Laird) [default], "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman) |


##### Output columns

| Column     | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
| `SNP`  | SNP identifier (e.g., rsID).                                                 |
| `CHR`  | Chromosome number where the SNP is located                                  |
| `BP`   | Base-pair position of the SNP on the chromosome                             |
| `P`    | P-value of the analysis/meta-analysis|
| `Z`    | Z-score value of the analysis/meta-analysis|
| `CI_lower`| Lower confidence interval|
| `CI_upper`| Upper confidence interval|
| `Weighted_Mean`|Analysis estimate|

## Plots 

You can generate Manhattan plots and QQ-plots from your GWAS analysis/meta-analysis results.

```bash
python plots.py   --i  ressults_file.txt
```

| Arguments                          | Description                                                                                 |
|-------------------------------|---------------------------------------------------------------------------------------------|
| `--i`, `--input`| Path to input tab separated file with SNP,CHR,BP,P columns       |
| `--o`, `--output`|  Prefix for output files (default: gwas_plot)                                                                    |
| `--dpi`         | Resolution for figures in dpi (default: 600)                                                        |


## Enrichment Analysis 

You can perform Manhattan Enrichment Analysis with gProfiler for your GWAS analysis/meta-analysis results.

```bash
python enrichment_analysis.py -i input_snps.txt -o gprofiler_results.csv -p gprofiler_manhattan.png -s 0.01 --organism hsapiens --sources GO:BP KEGG REAC --no_iea
```

| Arguments          | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `-i`, `--input`        | Input tab-delimited file with `'SNP'` and `'P'` columns                    |
| `-o`, `--output`       | Output CSV file for enrichment results                                     |
| `-p`, `--plot`         | Output PNG file for Enrichment Analysis Manhattan plot                                   |
| `-s`, `--significance` | Adjusted p-value significance threshold [default=0.05]                                   |
| `--organism`           | Organism name for g:Profiler (e.g., `'hsapiens'`)                          |
| `--sources`            | Data sources to query (e.g., GO:BP GO:CC KEGG). If omitted, all available sources will be used. |
| `--no_iea`             | Exclude electronic GO annotations (IEA)                                    |



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
