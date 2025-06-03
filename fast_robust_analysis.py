import os
import pandas as pd
import numpy as np
import math
import discrete_model, robust
import argparse
from scipy.stats import norm, cauchy
from decimal import Decimal, getcontext

# Set the precision for Decimal calculations
getcontext().prec = 1000


def cauchy_cdf(x, x0=0, gamma=1):
    """Cumulative distribution function (CDF) of the Cauchy distribution."""
    return 0.5 + (np.arctan((x - x0) / gamma) / np.pi)






def fast_robust_analysis(data,effect_size_type):
    
    data = data.dropna()
    # Combine dataframes into the first one
    file =  data
    file.index = file['SNP']
    snps = np.array(file.index)

    file = file.drop(['SNP'], axis=1)

    snp_keys = file.index.unique()
    snp_hash = {snp: file.loc[snp].values.tolist() for snp in snp_keys}

    p_value_min_p = []
    p_value_cauchy = []
    z_dom_list = []
    z_add_list = []
    z_rec_list = []
    
    p_dom = []
    p_rec = []
    p_add = []

    # Iterate through each SNP
    for snp_name in snp_keys:
        snp_data = pd.DataFrame(snp_hash[snp_name])
        #print(snp_name,snp_data)
        AA1, AB1, BB1, AA0, AB0, BB0 = [np.array(snp_data[0][col]) for col in range(2, 8)]

        aa1, ab1, bb1, aa0, ab0, bb0 = AA1, AB1, BB1, AA0, AB0, BB0
       
        if any([x == 0 for x in [aa1.astype(float), ab1.astype(float), bb1.astype(float), aa0.astype(float), ab0.astype(float), bb0.astype(float)]]):
            aa1 = aa1.astype(float) + 0.5
            ab1 = ab1.astype(float) + 0.5
            bb1 = bb1.astype(float) + 0.5
            aa0 = aa0.astype(float) + 0.5
            ab0 = ab0.astype(float) + 0.5
            bb0 = bb0.astype(float) + 0.5

        row_list = [aa1, ab1, bb1, aa0, ab0, bb0]
        #print(row_list)
        effect_dom, var_dom = discrete_model.model_dominant(row_list, effect_size_type)
        effect_add, var_add = discrete_model.model_additive(row_list, effect_size_type)
        effect_rec, var_rec = discrete_model.model_recessive(row_list, effect_size_type)

        z_dom = effect_dom / math.sqrt(var_dom) if var_dom > 0 else 0
        z_add = effect_add / math.sqrt(var_add) if var_add > 0 else 0
        z_rec = effect_rec / math.sqrt(var_rec) if var_rec > 0 else 0

        z_dom_list.append(z_dom)
        z_add_list.append(z_add)
        z_rec_list.append(z_rec)

        # Combine p-values using MinP and Cauchy
        p_vals = [norm.sf(abs(z)) * 2 for z in [z_dom, z_add, z_rec]]
        p_vals = np.array(p_vals)

        T = np.tan((0.5 - p_vals) * np.pi)
        t = np.sum(T) / 3

        # Calculate the combined p-value using the Cauchy distribution
        p_cauchy = Decimal(cauchy.sf(t))
        p_min = min(p_vals)
        combined_p = 1 - (1 - Decimal(p_min)) ** 3

        p_value_min_p.append(np.float64(combined_p))
        p_value_cauchy.append(np.float64(p_cauchy))
       
        p_dom.append(p_vals[0])
        p_rec.append(p_vals[2])
        p_add.append(p_vals[1])
        
        
    p_value_min_p = np.array(p_value_min_p)
    p_value_cauchy = np.array(p_value_cauchy)

    # Compute the tangent terms
    tan_terms = np.tan((0.5 - p_value_cauchy) * np.pi) + np.tan((0.5 - p_value_min_p) * np.pi)
    # Compute the CMC values
    p_cmc = 0.5 - np.arctan(tan_terms / 2) / np.pi

    # Compute the MCM values
    min_values = np.minimum(p_value_cauchy, p_value_min_p)
    p_mcm = np.minimum(1, 2 * min_values)

    return pd.DataFrame(
        {'SNP': snps,
         'CHR': np.array(file['CHR']),
         'BP': np.array(file['BP']),
         
         
         'Z_Dom': z_dom_list,
         'Z_Add': z_add_list,
         'Z_Rec': z_rec_list,
         
         'P_Dom': np.array(p_dom),
         'P_Add': np.array(p_add),
         'P_Rec': np.array(p_rec),
         
         'P_MinP': p_value_min_p,
         'P_CCT': p_value_cauchy,
         'P_CMC': p_cmc,
         'P_MCM': p_mcm

         })









# def fast_robust_analysis_sep(path, file_list, effect_size_type):
    # study_list = [pd.read_csv(os.path.join(path, file), sep='\t', encoding='latin1') for file in file_list]

    # # Combine dataframes into the first one
    # file = study_list[0]
    # file.index = file['SNP']
    # snps = np.array(file.index)

    # file = file.drop(['SNP'], axis=1)

    # snp_keys = file.index.unique()
    # snp_hash = {snp: file.loc[snp].values.tolist() for snp in snp_keys}

    # p_value_min_p = []
    # p_value_cauchy = []
    # z_dom_list = []
    # z_add_list = []
    # z_rec_list = []

    # # Iterate through each SNP
    # for snp_name in snp_keys:
        # snp_data = pd.DataFrame(snp_hash[snp_name])

        # AA1, AB1, BB1, AA0, AB0, BB0 = [np.array(snp_data[0][col]) for col in range(2, 8)]

        # aa1, ab1, bb1, aa0, ab0, bb0 = AA1, AB1, BB1, AA0, AB0, BB0

        # if any([x == 0 for x in [aa1, ab1, bb1, aa0, ab0, bb0]]):
            # aa1 = aa1.astype(float) + 0.5
            # ab1 = ab1.astype(float) + 0.5
            # bb1 = bb1.astype(float) + 0.5
            # aa0 = aa0.astype(float) + 0.5
            # ab0 = ab0.astype(float) + 0.5
            # bb0 = bb0.astype(float) + 0.5

        # row_list = [aa1, ab1, bb1, aa0, ab0, bb0]

        # effect_dom, var_dom = discrete_model.model_dominant(row_list, effect_size_type)
        # effect_add, var_add = discrete_model.model_additive(row_list, effect_size_type)
        # effect_rec, var_rec = discrete_model.model_recessive(row_list, effect_size_type)

        # z_dom = effect_dom / math.sqrt(var_dom) if var_dom > 0 else 0
        # z_add = effect_add / math.sqrt(var_add) if var_add > 0 else 0
        # z_rec = effect_rec / math.sqrt(var_rec) if var_rec > 0 else 0

        # z_dom_list.append(z_dom)
        # z_add_list.append(z_add)
        # z_rec_list.append(z_rec)

        # # Combine p-values using MinP and Cauchy
        # p_vals = [norm.sf(abs(z)) * 2 for z in [z_dom, z_add, z_rec]]
        # p_vals = np.array(p_vals)

        # T = np.tan((0.5 - p_vals) * np.pi)
        # t = np.sum(T) / 3

        # # Calculate the combined p-value using the Cauchy distribution
        # p_cauchy = Decimal(cauchy.sf(t))
        # p_min = min(p_vals)
        # combined_p = 1 - (1 - Decimal(p_min)) ** 3

        # p_value_min_p.append(np.float64(combined_p))
        # p_value_cauchy.append(np.float64(p_cauchy))
        
       
        
        
        
    # p_value_min_p = np.array(p_value_min_p)
    # p_value_cauchy = np.array(p_value_cauchy)

    # # Compute the tangent terms
    # tan_terms = np.tan((0.5 - p_value_cauchy) * np.pi) + np.tan((0.5 - p_value_min_p) * np.pi)
    # # Compute the CMC values
    # p_cmc = 0.5 - np.arctan(tan_terms / 2) / np.pi

    # # Compute the MCM values
    # min_values = np.minimum(p_value_cauchy, p_value_min_p)
    # p_mcm = np.minimum(1, 2 * min_values)

    # return pd.DataFrame(
        # {'SNP': snps,
         # 'CHR': np.array(file['CHR']),
         # 'BP': np.array(file['BP']),
         # 'Z_Dom': z_dom_list,
         # 'Z_Add': z_add_list,
         # 'Z_Rec': z_rec_list,
         
         # 'P_MinP': p_value_min_p,
         # 'P_CCT': p_value_cauchy,
         # 'P_CMC': p_cmc,
         # 'P_MCM': p_mcm
       
         # })


if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Fast robust analysis for genetic data.")

    parser.add_argument("--path", type=str, required=True, help="Path to the directory containing input files.")
    parser.add_argument("--file_list", type=str, required=True, nargs='+',
                        help="List of input data files (space-separated).")
    parser.add_argument("--effect_size_type", type=str, required=True, help="Type of effect size (e.g., 'OR', 'log').")
    parser.add_argument("--output", type=str, default="results.csv", help="Output file name for saving the results.")

    # Parse arguments
    args = parser.parse_args()

    # Perform the analysis
    final_results = fast_robust_analysis(path=args.path, file_list=args.file_list,
                                         effect_size_type=args.effect_size_type)

    # Save results to the output file
    final_results.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to {args.output}")
