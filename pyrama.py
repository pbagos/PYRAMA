import sys
import os
import argparse
import time
import pandas as pd
import polars as pl
import meta_analysis
import cont_meta_analysis
import fast_robust_analysis
import fast_robust_cont_analysis
import bayesian
import cont_bayesian
from typing import List, Optional


from concurrent.futures import ThreadPoolExecutor, as_completed


def process_snps(input_file: str, output_file: str) -> None:
    """
    Reads an input file, identifies SNPs occurring less than the maximum number of occurrences,
    and writes them to an output file.
    """
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    if 'SNP' not in df.columns:
        print("Error: The input file must contain a 'SNP' column.")
        return

    snp_counts = df['SNP'].value_counts()
    max_occurrences = snp_counts.max()
    filtered_snps = snp_counts[snp_counts < max_occurrences].index.tolist()

    try:
        with open(output_file, 'w') as f:
            for snp in filtered_snps:
                f.write(snp + '\n')
        print(f"Filtered SNPs have been written to {output_file}")
    except Exception as e:
        print(f"Error writing to output file: {e}")


def _read_tsv(file_path: str) -> pd.DataFrame:
    """
    Read a single TSV file into a DataFrame.
    Raises RuntimeError on failure, for centralized error handling.
    """
    try:
        df = pd.read_csv(file_path, sep='\t', dtype=str)
        return df
    except Exception as e:
        raise RuntimeError(f"Error reading file {file_path}: {e}")


 

def merge_input_files(
    file_list: List[str],
    max_workers: int = 1
) -> pd.DataFrame:
    """
    Reads and merges multiple TSV files by rows using Polars for fast I/O
    and an optional thread pool for file‐level concurrency.

    Parameters
    ----------
    file_list : List[str]
        Paths to TSV files.
    max_workers : int, optional
        Number of threads to use for reading files in parallel. Default is 1.

    Returns
    -------
    pandas.DataFrame
        Concatenated DataFrame sorted by 'SNP'.
    """
    # Helper to read a single TSV into a Polars DataFrame
    def _read_tsv_polars(fp: str) -> pl.DataFrame:
        return pl.read_csv(fp, separator='\t')

    polars_dfs: List[pl.DataFrame] = []
    # Parallel file reads
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_fp = {
            executor.submit(_read_tsv_polars, fp): fp
            for fp in file_list
        }
        for future in as_completed(future_to_fp):
            fp = future_to_fp[future]
            try:
                df_pl = future.result()
                polars_dfs.append(df_pl)
            except Exception as exc:
                print(f"Error processing {fp}: {exc}", file=sys.stderr)
                sys.exit(1)

    # Concatenate vertically and sort by 'SNP'
    merged_pl = pl.concat(polars_dfs, how='vertical').sort('SNP')

    # Convert back to pandas DataFrame
    merged_df = merged_pl.to_pandas()
    return merged_df


def gwas_meta_analysis(
    input_files: List[str],
    output_file: str,
    inheritance_model: Optional[str],
    effect_size_type: Optional[str],
    robust_method: Optional[str],
    type_of_effect: Optional[str],
    approximate_max: Optional[str],
    biv_ma: str = 'NO',
    imputation: bool = False,
    bayesian_meta: str = 'NO',
    r2threshold: Optional[str] = None,
    population: Optional[str] = None,
    maf: Optional[str] = None,
    ref: Optional[str] = None,
    imp_list: Optional[str] = None,
    nthreads: int = 1,
    missing_threshold: float = 0.5,
    het_est : str =' '
) -> None:
    """
    Run GWAS meta-analysis. Supports imputation when imputation=True.
    If missing_threshold is 0, imputation via pred_ld is performed without passing an imp_list file.
    """
    print(f"Running GWAS Meta-Analysis with {nthreads} threads")

    # Column definitions
    case_1_columns = ['SNP', 'CHR', 'BP', 'aa1', 'ab1', 'bb1', 'aa0', 'ab0', 'bb0']
    case_12_columns = ['SNP', 'CHR', 'BP', 'aa2', 'ab2', 'bb2', 'aa1', 'ab1', 'bb1', 'aa0', 'ab0', 'bb0']
    case_2_columns = ['SNP', 'CHR', 'BP', 'BETA', 'SE']
    case_22_columns = ['SNP', 'CHR', 'BP', 'BETA1', 'SE1', 'BETA2', 'SE2']
    case_3_columns = ['SNP', 'CHR', 'BP', 'xaa', 'sdaa', 'naa', 'xab', 'sdab', 'nab', 'xbb', 'sdbb', 'nbb']

    # Peek header of first file
    try:
        sample_header = pd.read_csv(input_files[0], sep='\t', dtype=str, nrows=0).columns.tolist()
    except Exception as e:
        print(f"Error reading header of {input_files[0]}: {e}")
        return

    # Case 2: direct BETA/SE
    if all(col in sample_header for col in case_2_columns):
        print("Detected BETA and SE in input files")
        if imputation:
            print(f"Missing threshold: {missing_threshold}")
            print("Performing imputation...")
            if not all([r2threshold, population, maf, ref]):
                raise ValueError("Imputation parameters are required when --imputation is enabled.")

            # Validate all files contain the correct case2 columns
            for file in input_files:
                header = pd.read_csv(file, sep='\t', dtype=str, nrows=0).columns.tolist()
                if not all(col in header for col in case_2_columns):
                    raise ValueError(f"File {file} does not match case_2_columns required for imputation.")

            # If missing_threshold is zero, skip computing imp_list files
            imputed_file_lists = []
            if float(missing_threshold) == 0:
                print("Missing threshold is 0: performing pred_ld without imp_list for all studies")
                for idx, file in enumerate(input_files):
                    output_dir = f"imputation_study{idx+1}"
                    os.makedirs(output_dir, exist_ok=True)

                    command = (
                        f"python3 pred_ld.py "
                        f"--file-path {file} "
                        f"--r2threshold {r2threshold} "
                        f"--pop {population} "
                        f"--maf {maf} "
                        f"--ref {ref} "
                    )
                    print(f"Running imputation for study {idx+1} with command:\n{command}\n")
                    os.system(command)
                    imputed_file_lists.append(file + "_imputation_results_chr_all.txt")

            else:
                # Read SNP lists and calculate missing SNPs
                dfs = [pd.read_csv(f, sep='\t', dtype=str)[['SNP']] for f in input_files]
                total_files = len(dfs)
                all_snps = pd.concat(dfs)
                snp_freqs = all_snps['SNP'].value_counts()
                min_required_studies = int(total_files * (1 - float(missing_threshold)))
                if min_required_studies == 0:
                    print("Cannot perform summary statistics imputation with this missing threshold. Lower the missing threshold")
                    return
                print(f"Minimum required studies for SNP presence: {min_required_studies} of {total_files}")

                rare_snps = snp_freqs[snp_freqs < min_required_studies].index.tolist()
                print(f"{len(rare_snps)} SNPs found below the threshold")

                imputed_file_lists = []
                # Generate missing lists and run pred_ld
                for idx, df in enumerate(dfs):
                    missing = list(set(rare_snps) - set(df['SNP'].unique()))
                    missing_file = f"missing_snps_study{idx+1}.txt"
                    with open(missing_file, 'w') as f:
                        for snp in missing:
                            f.write(snp + '\n')
                    print(f"Study {idx+1}: {len(missing)} missing SNPs written to {missing_file}")

                    output_dir = f"imputation_study{idx+1}"
                    os.makedirs(output_dir, exist_ok=True)

                    command = (
                        f"python3 pred_ld.py "
                        f"--file-path {input_files[idx]} "
                        f"--r2threshold {r2threshold} "
                        f"--pop {population} "
                        f"--maf {maf} "
                        f"--ref {ref} "
                        f"--imp_list {missing_file} "
                    )
                    print(f"Running imputation for study {idx+1} with command:\n{command}\n")
                    os.system(command)
                    
                    imputed_file_lists.append(input_files[idx] + "_imputation_results_chr_all.txt")

            # Meta-analysis on imputed files
            script_inputs = ' '.join(imputed_file_lists)
            start_time = time.time()
            
            os.system(f"./pyrama_beta_se_meta {nthreads} {het_est} {script_inputs} > {output_file}")
            elapsed = time.time() - start_time
            print(f"Elapsed time: {elapsed:.2f} seconds")
            print(f"Meta-analysis results saved to {output_file}")
         
            return
        
        # No imputation
        script_inputs = ' '.join(input_files)
        print("Running Meta-Analysis with BETA and SE...")
     
        start_time = time.time()
        os.system(f"./pyrama_beta_se_meta {nthreads} {het_est} {script_inputs} > {output_file}")
        elapsed = time.time() - start_time
        print(f"Elapsed time: {elapsed:.2f} seconds")
        print(f"Meta-analysis results saved to {output_file}")
        return

    # Merge for other cases
    merged_df = merge_input_files(input_files, max_workers=nthreads)
 
    # Bivariate allele counts
    if all(col in merged_df.columns for col in case_12_columns):
        print("Bivariate meta-analysis (allele counts)")
        data_subset = merged_df[case_12_columns]
        result = bivariate.biv_meta_analysis(
            data_subset, inheritance_model, effect_size_type,
            robust_method, type_of_effect, approximate_max
        )
        result.to_csv(output_file, sep='\t', index=False)

    # Bivariate GWAS BETA/SE
    elif all(col in merged_df.columns for col in case_22_columns):
        print("Bivariate meta-analysis with BETA and SE")
        data_subset = merged_df[case_22_columns]
        result = bivariate_gwas.beta_SE_meta(data_subset)
        result.to_csv(output_file, sep='\t', index=False)

    # Case 1: discrete counts
    elif all(col in merged_df.columns for col in case_1_columns):
        print("Discrete phenotype input detected")
        data_subset = merged_df[case_1_columns]
        if robust_method == 'FAST':
            print("Fast Robust methods analysis/meta-analysis")
            result = fast_robust_analysis.fast_robust_analysis(data_subset, effect_size_type, het_est)
        elif bayesian_meta == 'YES':
            print("Bayesian meta-analysis")
            result = bayesian.meta_analysis(
                data_subset, inheritance_model, effect_size_type,
                robust_method, approximate_max
            )
        else:
            print("Standard meta-analysis")
            result = meta_analysis.meta_analysis(
                data_subset, inheritance_model, effect_size_type,
                robust_method, type_of_effect, approximate_max
            )
        result.to_csv(output_file, sep='\t', index=False)

    # Case 3: continuous phenotype
    elif all(col in merged_df.columns for col in case_3_columns):
        print("Continuous phenotype input detected")
        data_subset = merged_df[case_3_columns]
        if bayesian_meta =='YES':
                print("Bayesian meta-analysis")
                result = cont_bayesian.meta_analysis(data_subset, inheritance_model, effect_size_type,robust_method, approximate_max )
         
                result.to_csv(output_file, sep='\t', index=False)
                

        elif robust_method == 'FAST':
            result = fast_robust_cont_analysis.fast_robust_analysis(data_subset , het_est)
        else:
            result = cont_meta_analysis.meta_analysis(data_subset, inheritance_model, robust_method, type_of_effect)
        result.to_csv(output_file, sep='\t', index=False)

    else:
        raise ValueError("Data does not match the required columns for any supported case.")

    print(f"Meta-analysis results saved to {output_file}")


if __name__ == "__main__":
    print(r"""
_______  __      __  _______    ______   __       __    ______  
|       \|  \    /  \|       \  /      \ |  \     /  \ /      \ 
| #######\\##\  /  ##| #######\|  ######\| ##\   /  ##|  ######\
| ##__/ ## \##\/  ## | ##__| ##| ##__| ##| ###\ /  ###| ##__| ##
| ##    ##  \##  ##  | ##    ##| ##    ##| ####\\  #### ##    ##
| #######    \####   | #######\| ########| ##\## ## ##| ########
| ##         | ##    | ##  | ##| ##  | ##| ## \###| ##| ##  | ##
| ##         | ##    | ##  | ##| ##  | ##| ##  \# | ##| ##  | ##
 \##          \##     \##   \## \##   \## \##      \## \##   \##
    """)
    print("PYRAMA: A Python tool for Robust Analysis and Meta-Analysis of genome-wide association studies")
    version = '1.0.0'
    print("Version " + version + "; July 2025")
    print("Copyright (C) Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("-------------------------------------------------------------------------------")

    parser = argparse.ArgumentParser(description="Perform GWAS Meta-Analysis.")
    parser.add_argument(
        "--i", "--input", "--input_file", nargs='+', required=True,
        help="Paths to input data files (space-separated)."
    )
    parser.add_argument(
        "--o", "--output", "--output_file", required=True,
        help="Path to save the output results."
    )
    parser.add_argument("--inheritance_model", help="Inheritance model: ADDITIVE, RECESSIVE or DOMINANT")
    parser.add_argument("--effect_size_type", help="Type of effect size: OR or CATT")
    parser.add_argument("--robust_method", help="Robust method: MIN, MAX or MERT, or FAST")
    parser.add_argument("--type_of_effect", help="Type of effect: FIXED or RANDOM")
    parser.add_argument("--approximate_max", help="Approximate maximum: YES or NO")
    parser.add_argument("--biv_ma", default="NO", help="Bivariate meta-analysis (default: NO)")
    parser.add_argument("--imputation", action="store_true", help="Enable imputation step")
    parser.add_argument("--bayesian_meta", default="NO", help="Bayesian meta-analysis (default: NO)")
    parser.add_argument("--het_est", required = False, default=" ", help="Heterogeneity estimator: 'DL' (DerSimonian-Laird) [default], "
                             "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman)")



    parser.add_argument("--r2threshold", help="R2 threshold for imputation.")
    parser.add_argument("--population", help="Population for imputation.")
    parser.add_argument("--maf", help="Minor allele frequency for imputation.")
    parser.add_argument("--ref", help="Reference panel for imputation.")
    parser.add_argument("--imp_list", help="List of SNPs to impute (optional).")
    parser.add_argument(
    "--missing_threshold", type=float, default=0.5,
    help="Missingness threshold (default: 0.5). E.g. 0.5 means SNP must exist in at least 50%  of the studies."
)

    parser.add_argument(
        "--nthreads", "-n", type=int, default=1,
        help="Number of parallel threads to use (default: 1)."
    )

    args = parser.parse_args()

    gwas_meta_analysis(
        args.i, args.o, args.inheritance_model, args.effect_size_type,
        args.robust_method, args.type_of_effect, args.approximate_max,
        args.biv_ma, args.imputation, args.bayesian_meta,
        args.r2threshold, args.population, args.maf, args.ref,
        args.imp_list, args.nthreads,args.missing_threshold,args.het_est

    )
