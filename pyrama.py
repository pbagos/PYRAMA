import meta_analysis
import beta_se_meta2
import sys
import os
import dask.dataframe as dd
import argparse
import imputation
import cont_meta_analysis
import fast_robust_analysis
import bayesian
import bivariate
import bivariate_gwas
import pandas as pd


def process_snps(input_file, output_file):
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


def merge_input_files(file_list):
    """
    Reads and merges multiple tab-delimited files by rows using pandas.
    """
    dataframes = []
    for file in file_list:
        try:
            df = pd.read_csv(file, sep='\t', dtype=str)
            dataframes.append(df)
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            sys.exit(1)

    merged_df = pd.concat(dataframes, axis=0, ignore_index=True)
    merged_df = merged_df.sort_values(by='SNP', ascending=True)

    
    return merged_df


def gwas_meta_analysis(input_files, output_file, inheritance_model, effect_size_type,
                       robust_method, type_of_effect, approximate_max, biv_ma='NO',
                       imputation=False, bayesian_meta = 'NO',r2threshold=None, population=None, maf=None, ref=None, imp_list=None):
    print("Running GWAS Meta-Analysis:")

    # Merge input files
    merged_df = merge_input_files(input_files)
    data_df =  merged_df
    print(data_df)
    case_1_columns = ['SNP', 'CHR', 'BP', 'aa1', 'ab1', 'bb1', 'aa0', 'ab0', 'bb0']
    case_12_columns = ['SNP', 'CHR', 'BP','aa2', 'ab2', 'bb2', 'aa1', 'ab1', 'bb1', 'aa0', 'ab0', 'bb0']
    case_2_columns = ['SNP', 'CHR', 'BP', 'BETA', 'SE']
    case_22_columns = ['SNP', 'CHR', 'BP', 'BETA1', 'SE1','BETA2', 'SE2']

    case_3_columns = ['SNP', 'CHR', 'BP', 'xaa', 'sdaa', 'naa', 'xab', 'sdab', 'nab', 'xbb', 'sdbb', 'nbb']
    
    if  all(col in data_df.columns for col in case_12_columns):
      
          data_subset = data_df[case_12_columns]  
          print("Bivariate meta-analysis")
          result = bivariate.biv_meta_analysis(data_subset,inheritance_model,effect_size_type,robust_method,type_of_effect,approximate_max)
          result.to_csv(output_file, sep='\t', index=False)
    
    
    elif  all(col in data_df.columns for col in case_22_columns):
          data_subset = data_df[case_22_columns]  
          print("Bivariate meta-analysis with BETA and SE")
          
          result = bivariate_gwas.beta_SE_meta(data_subset)
          
          result.to_csv(output_file, sep='\t', index=False)
          
          
    elif all(col in data_df.columns for col in case_1_columns):
        data_subset = data_df[case_1_columns]
           
        print("Discrete phenotype input detected")   
        if (robust_method =='FAST'):
            print("Fast Robust methods analysis/meta-analysis")
            result = fast_robust_analysis.fast_robust_analysis(data_subset,effect_size_type)
        if (bayesian_meta == 'YES'):
        
            print("Bayesian meta-analysis")
            result = bayesian.meta_analysis(data_subset, inheritance_model, effect_size_type, robust_method, approximate_max)
    
        
        
        else:
            print("Standard meta-analysis")

            result = meta_analysis.meta_analysis(data_subset, inheritance_model, effect_size_type,
                                             robust_method, type_of_effect, approximate_max)
        
        
        
        result.to_csv(output_file, sep='\t', index=False)
        
        
    
             
    elif all(col in data_df.columns for col in case_2_columns):
    
        if imputation:
            print("Performing imputation...")
            if not all([r2threshold, population, maf, ref]):
                raise ValueError("Imputation parameters are required when --imputation is enabled.")

            temp_input = "temp_merged_input.txt"
            merged_df.to_csv(temp_input, sep='\t', index=False)

            if imp_list:
                process_snps(temp_input, "snps_missing_from_studies.txt")
                imp_list = "snps_missing_from_studies.txt"
                os.system(f"python3 imputation.py {temp_input} results {r2threshold} {population} {maf} {ref} {imp_list}")
            else:
                os.system(f"python3 imputation.py {temp_input} results {r2threshold} {population} {maf} {ref}")

            input_file = "results/imputation_results.txt"
        else:
            input_file = "temp_merged_input.txt"
            merged_df.to_csv(input_file, sep='\t', index=False)

        print("Running Meta-Analysis with BETA and SE...")
        os.system(f"./pyrama_beta_SE_meta {input_file} > {output_file}")

    elif all(col in data_df.columns for col in case_3_columns):
        print("Continous phenotype input detected")
        data_subset = data_df[case_3_columns]
    
        result = cont_meta_analysis.meta_analysis(data_subset, inheritance_model, robust_method, type_of_effect)
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
    parser.add_argument("--i", "--input", "--input_file", nargs='+', required=True, help="Paths to input data files (space-separated).")
    parser.add_argument("--o", "--output", "--output_file", required=True, help="Path to save the output results.")
    parser.add_argument("--inheritance_model", required=False, help="Inheritance model to use. ADDITIVE, RECESSIVE or DOMINANT")
    parser.add_argument("--effect_size_type", required=False, help="Type of effect size. OR or CATT")
    parser.add_argument("--robust_method", required=False, help="Robust method to use. MIN, MAX or MERT, or FAST")
    parser.add_argument("--type_of_effect", required=False, help="Type of effect. FIXED or RANDOM")
    parser.add_argument("--approximate_max", required=False, help="Approximate maximum. YES or NO")
    parser.add_argument("--biv_ma", default="NO", help="Bivariate meta-analysis (default: NO).")
   


    parser.add_argument("--imputation", action="store_true", help="Enable imputation step (default: disabled).")
    parser.add_argument("--bayesian_meta",default="NO", help="Bayesian meta-analysis (default: NO).")
    parser.add_argument("--r2threshold", required=False, help="R2 threshold for imputation.")
    parser.add_argument("--population", required=False, help="Population for imputation.")
    parser.add_argument("--maf", required=False, help="Minor allele frequency for imputation.")
    parser.add_argument("--ref", required=False, help="Reference panel for imputation.")
    parser.add_argument("--imp_list", required=False, help="List of SNPs to impute (optional).")

    args = parser.parse_args()

    gwas_meta_analysis(args.i, args.o, args.inheritance_model, args.effect_size_type,
                       args.robust_method, args.type_of_effect, args.approximate_max, args.biv_ma,
                       args.imputation, args.bayesian_meta, args.r2threshold, args.population, args.maf, args.ref, args.imp_list)
