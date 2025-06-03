import meta_analysis
import beta_se_meta2
import sys
import os
import dask.dataframe as dd
import argparse
import imputation
import pandas as pd






def process_snps(input_file, output_file):
    """
    Reads an input file, identifies SNPs occurring less than the maximum number of occurrences,
    and writes them to an output file.

    Args:
        input_file (str): Path to the input file containing SNP data.
        output_file (str): Path to the output file where filtered SNPs will be saved.
    """
    # Read the input file into a DataFrame
    try:
        df = pd.read_csv(input_file, sep='\t')  # Adjust the delimiter if necessary
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # Ensure the input file has the required 'SNP' column
    if 'SNP' not in df.columns:
        print("Error: The input file must contain a 'SNP' column.")
        return

    # Count occurrences of each SNP
    snp_counts = df['SNP'].value_counts()

    # Find the maximum number of occurrences
    max_occurrences = snp_counts.max()

    # Filter SNPs with occurrences less than the maximum
    filtered_snps = snp_counts[snp_counts < max_occurrences].index.tolist()

    # Write the filtered SNPs to the output file
    try:
        with open(output_file, 'w') as f:
            for snp in filtered_snps:
                f.write(snp + '\n')
        print(f"Filtered SNPs have been written to {output_file}")
    except Exception as e:
        print(f"Error writing to output file: {e}")

# # Example usage
# input_file = "input.txt"  # Replace with your input file path
# output_file = "tab.txt"   # Replace with your desired output file path
# 
# process_snps(input_file, output_file)
# 




def gwas_meta_analysis(input_file, output_file, inheritance_model, effect_size_type,
                       robust_method, type_of_effect, approximate_max, biv_ma='NO',
                       imputation=False, r2threshold=None, population=None, maf=None, ref=None, imp_list=None):
    print("Running GWAS Meta-Analysis:")
    data_df = dd.read_csv(input_file, sep='\t', dtype={'ab0': 'float64', 'bb0': 'float64'})
    case_1_columns = ['SNP', 'CHR', 'BP', 'aa1', 'ab1', 'bb1', 'aa0', 'ab0', 'bb0']
    case_2_columns = ['SNP', 'CHR', 'BP', 'BETA', 'SE']
    case_3_columns = ['SNP', 'CHR', 'BP', 'xaa', 'sdaa', 'naa', 'xab', 'sdab', 'nab', 'xbb', 'sdbb', 'nbb']

    if all(col in data_df.columns for col in case_1_columns):
        data_subset = data_df[case_1_columns]
        data_subset = data_subset.compute()
        result = meta_analysis.meta_analysis(data_subset, inheritance_model, effect_size_type,
                                             robust_method, type_of_effect, approximate_max)
        result.to_csv(output_file, sep='\t', index=False)
    elif all(col in data_df.columns for col in case_2_columns):
        if imputation:
            print("Performing imputation...")
            if not all([r2threshold, population, maf, ref]):
                raise ValueError("Imputation parameters are required when --imputation is enabled.")
            
            if imp_list:  # Check if imp_list is provided
                process_snps(input_file, "snps_missing_from_studies.txt")
                imp_list = "snps_missing_from_studies.txt"
                os.system(f"python3 imputation.py {input_file} results {r2threshold} {population} {maf} {ref} {imp_list}")
            else:
                os.system(f"python3 imputation.py {input_file} results {r2threshold} {population} {maf} {ref}")
            
            input_file = "results/imputation_results.txt"  # Update the input_file to point to the imputed file
            
        print("Running Meta-Analysis with BETA and SE...")
        # Keep the columns from the input_file (SNP,CHR,BP,BETA,SE) and ensure that SNP is char,chr and bp are integers and BETA and SE are float

        os.system(f"./altmeta_fast {input_file} > {output_file}")
    elif all(col in data_df.columns for col in case_3_columns):
        data_subset = data_df[case_3_columns]
        data_subset = data_subset.compute()
        result = meta_analysis.meta_analysis(data_subset, inheritance_model, effect_size_type,
                                             robust_method, type_of_effect, approximate_max)
        result.to_csv(output_file, sep='\t', index=False)
    else:
        raise ValueError("Data does not match the required columns for Case 1, Case 2, or Case 3.")

    print(f"Meta-analysis results saved to {output_file}")

if __name__ == "__main__":
    print(r"""
_______  __      __  _______    ______   __       __   ______  
|       \|  \    /  \|       \  /      \ |  \     /  \ /      \ 
| #######\\##\  /  ##| #######\|  ######\| ##\   /  ##|  ######\
| ##__/ ## \##\/  ## | ##__| ##| ##__| ##| ###\ /  ###| ##__| ##
| ##    ##  \##  ##  | ##    ##| ##    ##| ####\  ####| ##    ##
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
    parser.add_argument("--i", "--input","--input_file", required=True, help="Path to the input data file.")
    parser.add_argument("--o", "--output", "--output_file",required=True, help="Path to save the output results.")
    parser.add_argument("--inheritance_model", required=False, help="Inheritance model to use. ADDITIVE, RECESSIVE or DOMINANT")
    parser.add_argument("--effect_size_type", required=False, help="Type of effect size. OR or CATT")
    parser.add_argument("--robust_method", required=False, help="Robust method to use. MIN, MAX or MERT")
    parser.add_argument("--type_of_effect", required=False, help="Type of effect. FIXED or RANDOM")
    parser.add_argument("--approximate_max", required=False, help="Approximate maximum. YES or NO")
    parser.add_argument("--biv_ma", default="NO", help="Bivariate meta-analysis (default: NO).")
    parser.add_argument("--imputation", action="store_true", help="Enable imputation step (default: disabled).")
    parser.add_argument("--r2threshold", required=False, help="R2 threshold for imputation.")
    parser.add_argument("--population", required=False, help="Population for imputation.")
    parser.add_argument("--maf", required=False, help="Minor allele frequency for imputation.")
    parser.add_argument("--ref", required=False, help="Reference panel for imputation.")
    parser.add_argument("--imp_list", required=False, help="List of SNPs to impute (default: disabled).")

    args = parser.parse_args()

    gwas_meta_analysis(args.i, args.o, args.inheritance_model, args.effect_size_type,
                       args.robust_method, args.type_of_effect, args.approximate_max, args.biv_ma,
                       args.imputation, args.r2threshold, args.population, args.maf, args.ref, args.imp_list)
