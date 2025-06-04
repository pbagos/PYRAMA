import argparse
import os
import pandas as pd
import numpy as np
import math
import robust,discrete_model
import sys
from scipy.integrate import quad
from scipy.stats import norm, chi2
from bivariate_gwas import mmom_multi
global new_weight,results_list,effect_size1,effect_size2,var1,var2

results_list = []



def biv_meta_analysis(data,inheritance_model,effect_size_type,robust_method,type_of_effect,approximate_max):
    
    file = data
    sumW = 0.0
    sumWY = 0.0
    sumWSquared = 0.0
    sumWYSquared = 0.0

    # get the unique SNPs
    # print(file)
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
    #          aa2  ab2   bb2   aa1  ab1   bb1  aa0  ab0   bb0
    # snp
    # rs111   40  587  1325   72  684  2180...
    # rs112   97  620  1100   88  579  1958...
    # rs113  110   97   845   69  224   845...
    # rs111  110   96  1230   63  777  1078..
    # rs112  963  441    67  126  897   903...
    #       aa1
    # 0     40
    # 1    110

    # Take each column (genotype)  For each SNP
    for snp_name in snp_keys:
        snp_name_list = []
        chrom_list = []

        pos_list = []

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
            AA2 = [np.array(snp_df[0][2])]
            AB2 = [np.array(snp_df[0][3])]
            BB2 = [np.array(snp_df[0][4])]
            AA1 = [np.array(snp_df[0][5])]
            AB1 = [np.array(snp_df[0][6])]
            BB1 = [np.array(snp_df[0][7])]
            AA0 = [np.array(snp_df[0][8])]
            AB0 = [np.array(snp_df[0][9])]
            BB0 = [np.array(snp_df[0][10])]





        else:
            num_of_studies = snp_df.shape[0]
            chrom = np.array(snp_df[0])
            pos = np.array(snp_df[1])
            AA2 = np.array(snp_df[2])
            AB2 = np.array(snp_df[3])
            BB2 = np.array(snp_df[4])
            AA1 = np.array(snp_df[5])
            AB1 = np.array(snp_df[6])
            BB1 = np.array(snp_df[7])
            AA0 = np.array(snp_df[8])
            AB0 = np.array(snp_df[9])
            BB0 = np.array(snp_df[10])

        # print("************")
        # print(AA1)
        # print(chrom)

        es_list = []
        var_list = []
        weight_list = []
        es_list2 = []
        var_list2 = []
        weight_list2 = []
        # for each study the SNP exists...
        for i in range(len(AA1)):

            chromosome = np.unique(chrom)[0]
            if num_of_studies == 1:
                aa2 = float(AA2[0])
                ab2 = float(AB2[0])
                bb2 = float(BB2[0])
                aa1 = float(AA1[0])
                ab1 = float(AB1[0])
                bb1 = float(BB1[0])
                aa0 = float(AA0[0])
                ab0 = float(AB0[0])
                bb0 = float(BB0[0])



            else:
                aa2 = float(AA2[i])
                ab2 = float(AB2[i])
                bb2 = float(BB2[i])
                aa1 = float(AA1[i])
                ab1 = float(AB1[i])
                bb1 = float(BB1[i])
                aa0 = float(AA0[i])
                ab0 = float(AB0[i])
                bb0 = float(BB0[i])

            # if any of the genotypes equals zero -> add 0.5 to every genotype

            if (aa2 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (ab2 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (bb2 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (aa1 == 0):


                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (ab1 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (bb1 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (aa0 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (ab0 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            if (bb0 == 0):
                aa2 += 0.5
                ab2 += 0.5
                bb2 += 0.5
                aa1 += 0.5
                ab1 += 0.5
                bb1 += 0.5
                aa0 += 0.5
                bb0 += 0.5
                ab0 += 0.5

            row_list = [aa2,ab2,bb2,aa1, ab1, bb1, aa0, ab0, bb0]
            effect_size = 0
            var = 0
            n2 = bb0 + bb1
            n1 = ab1 + ab0


            R1 = aa1 + ab1 + bb1
            R2 = aa2 + ab2 + bb2
            S = aa0 + ab0 + bb0

            N1 = R1 + S
            N2 = R2 + S

            total1 = N1
            total2 = N2
            pAA1 = (aa0 + aa1) / total1
            pAB1 = (ab0 + ab1) / total1
            pBB1 = (bb0 + bb1) / total1

            pAA2 = (aa0 + aa2) / total2
            pAB2 = (ab0 + ab2) / total2
            pBB2 = (bb0 + bb2) / total2


            corrDomAdd1 = robust.discreteCorrPrototypeFunction(pAA1, pAB1, pBB1)
            corrDomAdd2 = robust.discreteCorrPrototypeFunction(pAA2, pAB2, pBB2)

            # print (f'correcADD values:{snp_name,pBB, pAB, pAA}')

            corrRecAdd1 = robust.discreteCorrPrototypeFunction(pBB1, pAB1, pAA1)
            corrRecAdd2 = robust.discreteCorrPrototypeFunction(pBB2, pAB2, pAA2)

            corrRecDom1 = math.sqrt((pAA1 * pBB1) / ((1 - pAA1) * (1 - pBB1)))
            corrRecDom2 = math.sqrt((pAA2 * pBB2) / ((1 - pAA2) * (1 - pBB2)))

            row_list1 = [aa1, ab1, bb1, aa0, ab0, bb0]
            row_list2 = [aa2, ab2, bb2, aa0, ab0, bb0]

            if (inheritance_model == 'RECESSIVE'):
                effect_size1, var1 = discrete_model.model_recessive(row_list1, effect_size_type)
                effect_size2, var2 = discrete_model.model_recessive(row_list2, effect_size_type)

            # Discrete Dominant
            if (inheritance_model == 'DOMINANT'):
                # print("Discrete Dominant")

                effect_size1, var1 = discrete_model.model_dominant(row_list1, effect_size_type)
                effect_size2, var2 = discrete_model.model_dominant(row_list2, effect_size_type)


            # Discrete Additive
            if (inheritance_model == 'ADDITIVE'):
                effect_size1, var1 = discrete_model.model_additive(row_list1, effect_size_type)
                effect_size2, var2 = discrete_model.model_additive(row_list2, effect_size_type)

                # print('Discrete Additive')

            if (inheritance_model == 'ALL'):
                # Call the model functions
                effect_size_dom, var_dom = discrete_model.model_dominant(row_list, effect_size_type)
                effect_size_add, var_add = discrete_model.model_additive(row_list, effect_size_type)
                effect_size_rec, var_rec = discrete_model.model_recessive(row_list, effect_size_type)

                effect_size_dom2, var_dom2 = discrete_model.model_dominant(row_list2, effect_size_type)
                effect_size_add2, var_add2= discrete_model.model_additive(row_list2, effect_size_type)
                effect_size_rec2, var_rec2 = discrete_model.model_recessive(row_list2, effect_size_type)


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


                # print(f"Effect size and variance additive {effect_size_add,var_add}")
                if (math.sqrt(var_dom2)) == 0:
                    effectDom2 = 0
                else:
                    effectDom2 = effect_size_dom2 / math.sqrt(var_dom2)
                if (math.sqrt(var_rec2)) == 0:
                    effectRec2 = 0
                else:
                    effectRec2 = effect_size_rec2 / math.sqrt(var_rec2)

                if (math.sqrt(var_add2)) == 0:
                    effectAdd2 = 0
                else:
                    effectAdd2 = effect_size_add2 / math.sqrt(var_add2)


                # print(f'Effect additive {effectAdd}')
                weight = robust.robustWeight(R1, S, tSquared=0)
                weight2 = robust.robustWeight(R2, S, tSquared=0)

                # print(f"Weight with tSquared = 0  {weight}")
                if robust_method == 'MIN':
                    pChiSquared = robust.pValueChiSquared(row_list1)
                    effect_size = robust.DiscreteRobustMin(weight, effectAdd, pChiSquared)

                    pChiSquared2 = robust.pValueChiSquared(row_list2)
                    effect_size2 = robust.DiscreteRobustMin(weight2, effectAdd2, pChiSquared)

                if robust_method == 'MAX':
                    effect_size = robust.DiscreteRobustMax(effectDom, effectRec, effectAdd, corrDomAdd1, corrRecAdd1,
                                                           corrRecDom1,
                                                           weight, approximate_max)
                    effect_size2 = robust.DiscreteRobustMax(effectDom2, effectRec2, effectAdd2, corrDomAdd2, corrRecAdd2,
                                                           corrRecDom2,
                                                           weight2, approximate_max)

                if robust_method == 'MERT':
                    # print(robust_method)
                    effect_size = robust.DiscteteRobustMert(effectDom, effectRec, corrRecDom1, weight)
                    effect_size2 = robust.DiscteteRobustMert(effectDom2, effectRec2, corrRecDom2, weight2)

            # if inheritance_model != 'all':
            #     weight = 1 / var

            # get effect_sizes and variances from the model
            es_list.append(effect_size1)
            var_list.append(var1)
            es_list2.append(effect_size2)
            var_list2.append(var2)
            snp_name_list.append(snp_name)
            chrom_list.append(chromosome)
            pos_list.append(pos[0])
          
            
                # Combine ys1 and ys2 into a 2D array for input into mmom_multi
        ys = np.column_stack((es_list, es_list2))  # This creates a 2D array: each row [ys1, ys2]

        # Combine vars1 and vars2 into a 2D array (optional, depends on how your function uses them)
        vars_ = np.column_stack((var_list, var_list2))  # Create a 2D array for variances
    
        res = mmom_multi(ys, vars_)
        beta_hat = res['beta_hat']
        beta_cov = res['beta_cov']
        
        
        
        # Compute the Wald statistic and associated p-value.
        try:
            inv_cov = np.linalg.pinv(beta_cov)
            # Wald statistic: beta_hat' * inv_cov * beta_hat.
            wald = beta_hat.T @ inv_cov @ beta_hat
            p_value = chi2.sf(wald, df=beta_cov.shape[0])
        except np.linalg.LinAlgError:
            wald = np.nan
            p_value = np.nan
            
            
        results_list.append([snp_name,chromosome,pos[0],beta_hat[0],beta_hat[1],np.sqrt(np.diag(beta_cov))[0],np.sqrt(np.diag(beta_cov))[1],wald,p_value ])
        
    #print(pd.DataFrame(snp_name_list, chrom_list,pos_list,results_list,columns = ['SNP','CHR','BP','Beta_average','Beta_Standard_Error','p_value_beta','Delta_hat','Delta_Standard_Error','p_value_delta']))
    # print (results_list)
    #df1 = pd.DataFrame({'SNP':snp_name_list, 'CHR':chrom_list, 'BP':pos_list})
    # print(df1)
    df2 = pd.DataFrame(results_list,columns = ['SNP','CHR','BP','meta_BETA1','meta_BETA2','meta_SE1','meta_SE2','Wald' ,'P' ])
    #total = pd.concat([df1, df2], axis=1)
   
    return df2
  
  
  
  
# # Main execution
# if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Run bivariate meta-analysis for GWAS data.")
    # parser.add_argument("--path", required=True, help="Path to the directory containing GWAS files.")
    # parser.add_argument("--file_list", nargs='+', required=True, help="List of GWAS data files.")
    # parser.add_argument("--inheritance_model", required=True, choices=['RECESSIVE', 'DOMINANT', 'ADDITIVE', 'ALL'],
                        # help="Inheritance model to use.")
    # parser.add_argument("--effect_size_type", required=True, help="Type of effect size calculation.")
    # parser.add_argument("--robust_method", required=True, choices=['MIN', 'MAX', 'MERT'],
                        # help="Robust method to apply.")
    # parser.add_argument("--type_of_effect", required=True, help="Type of effect to compute.")
    # parser.add_argument("--approximate_max", required=True, help="Approximation max value.")
    # parser.add_argument("--output", required=True, help="Path and name of the output file (tab-separated).")

    # args = parser.parse_args()

    # # Run the function with command-line arguments
    # results = biv_meta_analysis(
        # path=args.path,
        # file_list=args.file_list,
        # inheritance_model=args.inheritance_model,
        # effect_size_type=args.effect_size_type,
        # robust_method=args.robust_method,
        # type_of_effect=args.type_of_effect,
        # approximate_max=args.approximate_max
    # )
    # output_file = args.output
    
    # # Save the results
    # results.to_csv(output_file, sep='\t', index=False)
    # print(f"Bivariate Meta-analysis completed. Results saved to {output_file}")
