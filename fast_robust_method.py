import os
import pandas as pd
import numpy as np
import math
import robust,discrete_model
import sys
import argparse
from decimal import Decimal, getcontext
from scipy.stats import cauchy
from scipy.integrate import quad
from scipy.stats import norm, chi2
from bayes_factor import bayes_factor

global new_weight,results_list,bf
# Set precision for Decimal calculations
getcontext().prec = 400  # Set a high precision


results_list = []



# Define a function to calculate p-value from MCM
def calculate_mcm(row):
    p_cau = row["P_CCT"]
    p_min = row["P_MinP"] 
    return min(1, 2 * min(p_cau, p_min))

# Define a function to calculate p-value from CMC
def calculate_cmc(row):
    p_cau = row["P_CCT"]
    p_min = row["P_MinP"] 
    return 0.5 - math.atan(
        sum([math.tan((0.5 - p) * math.pi) for p in [p_cau, p_min]]) / 2
    ) / math.pi

def cauchy_cdf(x, x0=0, gamma=1):
    """Cumulative distribution function (CDF) of the Cauchy distribution."""
    return 0.5 + (np.arctan((x - x0) / gamma) / np.pi)


def CCT(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values, dtype=np.float64)
    # Ensure that all p-values are within the valid range (0, 1)
    # p_values = np.clip(p_values, 1e-16, 1 - 1e-16)

    k = len(p_values)

    T = np.tan((0.5 - p_values) * np.pi)

    t = np.sum(T) / k
    # Calculate the combined p-value using the Cauchy distribution
    combined_p = 1 - cauchy_cdf(t)
    #combined_p = Decimal(cauchy.sf(t)*2)

    return combined_p


# MinP combination
def minP(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values, dtype=np.float64)
    
    
   
    # Ensure that all p-values are within the valid range (0, 1)
    # p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
    k = len(p_values)
 
    
    # Calculate the combined p-value as the minimum p-value
 
    min_p = Decimal(np.min(p_values))
    combined_p = Decimal(1) - (Decimal(1) - min_p) ** 3

    return np.float64(combined_p)

# Wrapper function to run meta-analysis for each model and collect p-values
def fast_robust(path, file_list, effect_size_type, type_of_effect):
    models = ['DOMINANT', 'RECESSIVE', 'ADDITIVE']
    model_pvalues = {model: {} for model in models}  # Dictionary to store p-values by model

    # Iterate over models and run meta-analysis
    for model in models:
        print(f"Running meta-analysis for inheritance model: {model}")
        results = meta_analysis(
            path=path,
            file_list=file_list,
            inheritance_model=model,
            effect_size_type=effect_size_type,
            type_of_effect=type_of_effect 
    
        )
        for index, row in results.iterrows():
            snp = row['SNP']
            p_value = row['P_Value']
            model_pvalues[model][snp] = p_value

    # Combine results into a single DataFrame
    snp_list = list(model_pvalues['DOMINANT'].keys())  # Assume all models share the same SNPs
    final_results = pd.DataFrame({
        'SNP': snp_list,
        'P_ADD': [model_pvalues['ADDITIVE'].get(snp, np.nan) for snp in snp_list],
        'P_DOM': [model_pvalues['DOMINANT'].get(snp, np.nan) for snp in snp_list],
        'P_REC': [model_pvalues['RECESSIVE'].get(snp, np.nan) for snp in snp_list],
    })

    # Extract chromosome (CHR) and base-pair position (BP) information from one of the studies
    study_list = [pd.read_csv(path + file, sep='\t', encoding='latin1') for file in file_list]
    reference_study = study_list[0]
    chr_bp_mapping = reference_study.set_index('SNP')[['CHR', 'BP']].to_dict()

    final_results['CHR'] = final_results['SNP'].map(chr_bp_mapping['CHR'])
    final_results['BP'] = final_results['SNP'].map(chr_bp_mapping['BP'])
    
    #MinP and CCT
    final_results['P_MinP'] = final_results[["P_ADD", "P_DOM", "P_REC"]].apply(minP, axis = 1)
    final_results['P_CCT'] = final_results[["P_ADD", "P_DOM", "P_REC"]].apply(CCT, axis = 1)
    
  
    
    # Apply the functions to calculate P_MCM and P_CMC
    final_results['P_MCM'] = final_results.apply(calculate_mcm, axis=1)
    final_results['P_CMC'] = final_results.apply(calculate_cmc, axis=1)



    # Rearrange columns
    final_results = final_results[['SNP', 'CHR', 'BP', 'P_ADD', 'P_DOM', 'P_REC','P_MinP','P_CCT',"P_MCM","P_CMC"]]
    return final_results

# Meta-analysis function for a single model
def meta_analysis(path, file_list, inheritance_model, effect_size_type, type_of_effect):
    results_list = []

    study_list = [pd.read_csv(path + file, sep='\t', encoding='latin1') for file in file_list]
    snp_keys = study_list[0]['SNP'].unique()

    for snp_name in snp_keys:
        sumWY = 0
        sumW = 0
        weight_list = []
        es_list = []

        for study in study_list:
            if snp_name in study['SNP'].values:
                snp_data = study[study['SNP'] == snp_name].iloc[0]
                row_list = [snp_data[f'aa1'], snp_data[f'ab1'], snp_data[f'bb1'],
                            snp_data[f'aa0'], snp_data[f'ab0'], snp_data[f'bb0']]
                row_list = [x + 0.5 if x == 0 else x for x in row_list]

                if inheritance_model == 'DOMINANT':
                    effect_size, var = discrete_model.model_dominant(row_list, effect_size_type)
                elif inheritance_model == 'RECESSIVE':
                    effect_size, var = discrete_model.model_recessive(row_list, effect_size_type)
                elif inheritance_model == 'ADDITIVE':
                    effect_size, var = discrete_model.model_additive(row_list, effect_size_type)

                es_list.append(effect_size)
                weight = 1 / var if var != 0 else 0
                weight_list.append(weight)

                sumWY += effect_size * weight
                sumW += weight

        if sumW > 0:
            weighted_mean = sumWY / sumW
            z = abs(weighted_mean / math.sqrt(1 / sumW))
            p_value = norm.sf(abs(z)) * 2
        else:
            p_value = float('nan')

        results_list.append({'SNP': snp_name, 'P_Value': p_value})

    return pd.DataFrame(results_list)

if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Meta-analysis tool for genetic data.")
    
    parser.add_argument("--path", type=str, required=True, help="Path to the directory containing the data files.")
    parser.add_argument("--file_list", type=str, required=True, nargs='+', help="List of input data files.")
    parser.add_argument("--effect_size_type", type=str, required=True, help="Type of effect size (e.g., 'log', 'odds').")
    parser.add_argument("--type_of_effect", type=str, required=True, help="Type of effect (e.g., 'binary', 'continuous').")
    parser.add_argument("--output", type=str, default="results.csv", help="Output file name for the results.")

    # Parse arguments
    args = parser.parse_args()

    # Run the meta-analysis
    final_results = fast_robust(
        path=args.path,
        file_list=args.file_list,
        effect_size_type=args.effect_size_type,
        type_of_effect=args.type_of_effect
    )

    # Save results to the output file
    final_results.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to {args.output}")
# #Example
# df = fast_robust(path ='~/shiny-server/pyRAMA/', file_list = ['popi_m1.txt'],type_of_effect = 'RANDOM', effect_size_type = 'OR') 
# 
# # Convert entire data frame as string and print
# print(df.to_string())
