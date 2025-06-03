import argparse

import pandas as pd
import numpy as np
import math
import robust, discrete_model
from scipy.stats import norm, chi2

bayesian_results_list = []

def bayesian_meta_analysis_discrete(beta_list, se_list,n_i):
    y_i = beta_list
    ste= se_list
    z_list = []
    E_m_list = []
    V_t_list = []
    E_t_list = []
    V_m_list = []
    CI_low = []
    CI_up = []
    g_list = []
    p_values_list = []



    mean_y_i = np.mean(y_i)
    k = len(y_i)
    if k < 3:
        a = 2
        b = 2
    else:
        a = 0
        b = 2

    RSSb_first = 0
    for i in range(len(y_i)):
        RSSb_first += y_i[i] ** 2

    RSSb = RSSb_first - k * mean_y_i ** 2

    first_big_sum_q = []
    for i in range(len(y_i)):
        sens_part1 = (n_i[i] - 3)
        if sens_part1 == 0.0:
            first_big_sum_q.append(0)
        else:
            q1 = (n_i[i] * ste[i] * ste[i]) / sens_part1
            q2 = (mean_y_i * (k - 3) + y_i[i]) / k
            q3 = ((mean_y_i * mean_y_i * (mean_y_i - y_i[i]) + mean_y_i * y_i[i] * y_i[i]) * (
                        k + 2 * a + 1) * b) / (2 * (1 + b * RSSb / 2))
            q4 = q1 * (q2 - q3)
            first_big_sum_q.append(q4)
    first_big_sum = sum(first_big_sum_q)

    second_big_sum_q = []
    for i in range(len(y_i)):
        sens_part2 = (n_i[i] - 3)
        if sens_part2 == 0.0:
            second_big_sum_q.append(0)
        else:
            second_big_sum_q.append(((n_i[i] * ste[i] * ste[i]) / sens_part2) * (((k - 1) / k) - (
                        ((mean_y_i * (mean_y_i - y_i[i]) + y_i[i] ** 2) * (k + 2 * a + 1) * b) / (
                            2 * (1 + b * RSSb / 2)))))

    t1 = b * (k + 2 * a - 1)
    t2 = 2 * (1 + b * RSSb / 2)
    f = t1 / t2
    upper_part = mean_y_i - (f * first_big_sum)
    second_big_sum = sum(second_big_sum_q)
    lower_part = 1 - (f * second_big_sum)
    E_m = upper_part / lower_part
    sens_part2 = b * k * (2 * a + k - 3)
    if sens_part2 == 0:
        sens_part2 = 0.0000000000000000000000001
    sens_part3 = b * (2 * a + k - 3)
    if sens_part3 == 0:
        sens_part3 = 0.0000000000000000000000001
    V_mu = (2 * (1 + b * RSSb / 2)) / sens_part2
    E_tau_square = (2 * (1 + b * RSSb / 2)) / sens_part3
    sens_part4 = b ** 2 * (2 * a + k - 3) ** 2 * (k + 2 * a - 5)
    if sens_part4 == 0:
        sens_part4 = 0.000000001
    V_tau_square = (8 * (1 + b * RSSb / 2) ** 2) / sens_part4
    conf_int_up = E_m + 1.96 * math.sqrt(V_mu)
    conf_int_lower = E_m - 1.96 * math.sqrt(V_mu)
    z = E_m / math.sqrt(V_mu)
    p_value = norm.sf(abs(z)) * 2


    return [p_value,E_m,E_tau_square,V_mu,V_tau_square,conf_int_lower,conf_int_up,z]



def meta_analysis(path, file_list, inheritance_model, effect_size_type, robust_method, approximate_max):
    study_list = []
    for file in file_list:
        study_list.append(pd.read_csv(path + file, sep='\t', encoding='latin1'))
    # merge... files
    file = study_list[0]
    sumW = 0.0
    sumWY = 0.0
    sumWSquared = 0.0
    sumWYSquared = 0.0

    # get the unique SNPs

    snp_keys = file['SNP'].unique()
    # print(snp_keys)
    snp_values = []

    file.index = file['SNP']
    file = file.drop(['SNP'], axis=1)

    # print(file)

    snp_hash = {}

    for snp in snp_keys:
        snp_values.append(file.loc[snp].values.tolist())

    for i, snp in enumerate(snp_keys):
        snp_hash.update({snp: snp_values[i]})

    # print(snp_hash)

    # Inheritance models

    # For each SNP calculate the Odds Ratios

    # Discrete recessive

    # Odds ratio
    #             aa1  ab1   bb1  aa0  ab0   bb0
    # snp
    # rs111   40  587  1325   72  684  2180
    # rs112   97  620  1100   88  579  1958
    # rs113  110   97   845   69  224   845
    # rs111  110   96  1230   63  777  1078
    # rs112  963  441    67  126  897   903
    #       aa1
    # 0     40
    # 1    110

    # Take each column (genotype)  For each SNP
    for snp_name in snp_keys:
        sumWY = 0
        sumW = 0
        sumWSquared = 0
        sumWYSquared = 0
        snp_df = pd.DataFrame(snp_hash[snp_name])
        # print(snp_df.shape[1])
        if snp_df.shape[1] == 1:
            num_of_studies = 1
            chrom = [np.array(snp_df[0][0])]
            pos = [np.array(snp_df[0][1])]
            AA1 = [np.array(snp_df[0][2])]
            AB1 = [np.array(snp_df[0][3])]
            BB1 = [np.array(snp_df[0][4])]
            AA0 = [np.array(snp_df[0][5])]
            AB0 = [np.array(snp_df[0][6])]
            BB0 = [np.array(snp_df[0][7])]
        else:
            num_of_studies = snp_df.shape[0]
            chrom = np.array(snp_df[0])
            pos = np.array(snp_df[1])
            AA1 = np.array(snp_df[2])
            AB1 = np.array(snp_df[3])
            BB1 = np.array(snp_df[4])
            AA0 = np.array(snp_df[5])
            AB0 = np.array(snp_df[6])
            BB0 = np.array(snp_df[7])

        # print("************")
        # print(AA1)
        # print(chrom)

        es_list = []
        var_list = []
        weight_list = []
        n_i_list = []
        # for each study the SNP exists...
        for i in range(len(AA1)):

            chromosome = np.unique(chrom)[0]
            if num_of_studies == 1:
                aa1 = float(AA1[0])
                ab1 = float(AB1[0])
                bb1 = float(BB1[0])
                aa0 = float(AA0[0])
                ab0 = float(AB0[0])
                bb0 = float(BB0[0])
            else:
                aa1 = float(AA1[i])
                ab1 = float(AB1[i])
                bb1 = float(BB1[i])
                aa0 = float(AA0[i])
                ab0 = float(AB0[i])
                bb0 = float(BB0[i])

            # Check if any of the genotypes equals zero
            if aa1 == 0 or ab1 == 0 or bb1 == 0 or aa0 == 0 or bb0 == 0 or ab0 == 0:
                # Add 0.5 to every genotype
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            row_list = [aa1, ab1, bb1, aa0, ab0, bb0]
            effect_size = 0
            var = 0
            n2 = bb0 + bb1
            n1 = ab1 + ab0

            R = aa1 + ab1 + bb1
            S = aa0 + ab0 + bb0
            N = R + S

            total = N
            pAA = (aa0 + aa1) / total
            pAB = (ab0 + ab1) / total
            pBB = (bb0 + bb1) / total
            corrDomAdd = robust.discreteCorrPrototypeFunction(pAA, pAB, pBB)
            # print (f'correcADD values:{snp_name,pBB, pAB, pAA}')

            corrRecAdd = robust.discreteCorrPrototypeFunction(pBB, pAB, pAA)

            corrRecDom = math.sqrt((pAA * pBB) / ((1 - pAA) * (1 - pBB)))

            if (inheritance_model == 'RECESSIVE'):
                effect_size, var = discrete_model.model_recessive(row_list, effect_size_type)

            # Discrete Dominant
            if (inheritance_model == 'DOMINANT'):
                # print("Discrete Dominant")

                effect_size, var = discrete_model.model_dominant(row_list, effect_size_type)

            # Discrete Additive
            if (inheritance_model == 'ADDITIVE'):
                effect_size, var = discrete_model.model_additive(row_list, effect_size_type)
                # print('Discrete Additive')

            if (inheritance_model == 'ALL'):
                # Call the model functions
                effect_size_dom, var_dom = discrete_model.model_dominant(row_list, effect_size_type)
                effect_size_add, var_add = discrete_model.model_additive(row_list, effect_size_type)
                effect_size_rec, var_rec = discrete_model.model_recessive(row_list, effect_size_type)
                # print(f"Effect size and variance additive {effect_size_add,var_add}")
                if (math.sqrt(var_dom)) == 0:
                    effectDom = 0
                else:
                    effectDom = effect_size_dom / math.sqrt(var_dom)
                if (math.sqrt(var_rec)) == 0:
                    effectRec = 0
                else:
                    effectRec = effect_size_rec / math.sqrt(var_rec)

                if (math.sqrt(var_add)) == 0:
                    effectAdd = 0
                else:
                    effectAdd = effect_size_add / math.sqrt(var_add)
                # print(f'Effect additive {effectAdd}')
                weight = robust.robustWeight(R, S, tSquared=0)
                # print(f"Weight with tSquared = 0  {weight}")
                if robust_method == 'MIN':
                    pChiSquared = robust.pValueChiSquared(row_list)
                    effect_size = robust.DiscreteRobustMin(weight, effectAdd, pChiSquared)

                if robust_method == 'MAX':
                    effect_size = robust.DiscreteRobustMax(effectDom, effectRec, effectAdd, corrDomAdd, corrRecAdd,
                                                           corrRecDom,
                                                           weight, approximate_max)

                if robust_method == 'MERT':
                    # print(robust_method)
                    effect_size = robust.DiscteteRobustMert(effectDom, effectRec, corrRecDom, weight)

            # if inheritance_model != 'all':
            #     weight = 1 / var

            # get effect_sizes and variances from the model
            es_list.append(effect_size)
            var_list.append(var)
            n_i_list.append((R * S) / N)

        bayesian_result = bayesian_meta_analysis_discrete(es_list,var_list,n_i_list)

        bayesian_result.insert(0, snp_name)# insert snp
        bayesian_result.insert(1, chromosome) # insert chrom
        bayesian_result.insert(2, int(pos[0]) ) # insert bp

 

        bayesian_results_list.append( bayesian_result)
    return pd.DataFrame(bayesian_results_list, columns=['SNP', 'CHR', 'BP', 'P',"E_m","E_tau_square","V_mu","V_tau_square","conf_int_lower","conf_int_up",'Z'])



def main():
  
    parser = argparse.ArgumentParser(description="Perform Bayesian meta-analysis on GWAS data.")
    parser.add_argument("--path", required=True, help="Path to the directory containing input files.")
    parser.add_argument("--file_list", required=True, nargs='+', help="List of input files separated by spaces.")
    parser.add_argument("--inheritance_model", required=True, choices=["RECESSIVE", "DOMINANT", "ADDITIVE", "ALL"],
                        help="Inheritance model to use.")
    parser.add_argument("--effect_size_type", required=True, help="Type of effect size to calculate.")
    parser.add_argument("--robust_method", required=True, choices=["MIN", "MAX", "MERT"],
                        help="Robust method to use for calculations.")
    parser.add_argument("--approximate_max", required=True,
                        help="Approximate maximum value for robust calculations.")
    parser.add_argument("--output", required=True, help="Path and name of the output file (tab-separated).")

    args = parser.parse_args()

    # Parse arguments
    path = args.path
    file_list = args.file_list
    inheritance_model = args.inheritance_model
    effect_size_type = args.effect_size_type
    robust_method = args.robust_method
    approximate_max = args.approximate_max
    output_file = args.output

    # Perform meta-analysis
    results = meta_analysis(path, file_list, inheritance_model, effect_size_type, robust_method, approximate_max)

    # Save results as a tab-separated text file
    results.to_csv(output_file, sep="\t", index=False)
    print(f"Bayesian Meta-analysis completed. Results saved to {output_file}")

if __name__ == "__main__":
    main()
