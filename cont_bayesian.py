import argparse

import pandas as pd
import numpy as np
import math
import cont_robust, cont_model
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
    if k < 3: # if studies are below 3, a=2, b=2
        a = 2
        b = 2
    else: # if studies are above 3, a=0, b=2
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



def meta_analysis(data, inheritance_model, effect_size_type, robust_method, approximate_max):


    file = data
 
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

    for snp_name in snp_keys:
            # print(snp_name)
            sumWY = 0
            sumW = 0
            sumWSquared = 0
            sumWYSquared = 0
            snp_df = pd.DataFrame(snp_hash[snp_name])
            # print(snp_df.shape[1])
            # print(f'Shape of rows : {len(snp_df[0])}')
             
            #print(snp_df)
            # If the study is single
            if snp_df.shape[1] == 1:

                num_of_studies = 1
                snp_df_single  = np.array(snp_df)  
                 
                chrom = np.array([snp_df[0][0]],dtype = int)
                pos = np.array([snp_df[0][1]],dtype = int)

                XAA = np.array([snp_df[0][2]],dtype = float)
                SDAA = np.array([snp_df[0][3]],dtype = float)
                NAA = np.array([snp_df[0][4]],dtype = float)

                XAB = np.array([snp_df[0][5]],dtype = float)
                SDAB = np.array([snp_df[0][6]],dtype = float)
                NAB = np.array([snp_df[0][7]],dtype = float)

                XBB = np.array([snp_df[0][8]],dtype = float)
                SDBB = np.array([snp_df[0][9]],dtype = float)
                NBB = np.array([snp_df[0][10]],dtype = float)

             

            # Multiple Studies
            else:


                num_of_studies = snp_df.shape[0]

                chrom = np.array(snp_df[0],dtype = int)
                pos = np.array(snp_df[1],dtype = int)

                XAA = np.array(snp_df[2],dtype = float)
                SDAA = np.array(snp_df[3],dtype = float)
                NAA = np.array(snp_df[4],dtype = float)

                XAB = np.array(snp_df[5],dtype = float)
                SDAB = np.array(snp_df[6],dtype = float)
                NAB = np.array(snp_df[7],dtype = float)

                XBB = np.array(snp_df[8],dtype = float)
                SDBB = np.array(snp_df[9],dtype = float)
                NBB = np.array(snp_df[10],dtype = float)

                # print('------------------')

                # print(XAA,SDAA,NAA)

            es_list = []
            var_list = []
            n_i_list = [] 
            weight_list = []
            # for each study the SNP exists...
         
            for i in range(len(XAA)):
                
                if num_of_studies == 1:
                    xaa = XAA [0]
                    sdaa = SDAA [0]
                    naa = NAA [0]

                    xab = XAB [0]
                    sdab = SDBB [0]
                    nab = NAB [0]

                    xbb =  XBB [0]
                    sdbb = SDBB [0]
                    nbb =  NBB [0]

                else:

                    xaa = XAA [i]
                    sdaa = SDAA [i]
                    naa = NAA [i]

                    xab = XAB [i]
                    sdab = SDBB [i]
                    nab = NAB [i]

                    xbb =  XBB [i]
                    sdbb = SDBB [i]
                    nbb =  NBB [i]

                row_list = [xaa,sdaa,naa,xab,sdab,nab,xbb,sdbb,nbb]
                # print('----------row_list----------')
                

                effect_size = 0
                var = 0

              
          
        
                N = naa +nab +nbb
                total = N
                # pAA = (aa0 + aa1) / total
                # pAB = (ab0 + ab1) / total
                # pBB = (bb0 + bb1) / total
                
                
              

                if (inheritance_model == 'RECESSIVE'):
                    effect_size, var = cont_model.cont_model_recessive(row_list)

                # Discrete Dominant
                if (inheritance_model == 'DOMINANT'):
                    # print("Discrete Dominant")

                    effect_size, var = cont_model.cont_model_dominant(row_list)

                # Discrete Additive
                if (inheritance_model == 'ADDITIVE'):
                    effect_size, var = cont_model.cont_model_additive(row_list )
                    # print('Discrete Additive')

                if (inheritance_model == 'ALL'):
                    # Call the model functions
                    effect_size_dom, var_dom = cont_model.cont_model_dominant(row_list )
                    effect_size_add, var_add = cont_model.cont_model_additive(row_list )
                    effect_size_rec, var_rec = cont_model.cont_model_recessive(row_list )
                   

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
                    weight = cont_robust.robustWeight(row_list, tSquared=0)
                    # print(f"Weight with tSquared = 0  {weight}")
            
                    if robust_method == 'MAX':
                        effect_size = cont_robust.ContRobustMax(effectDom, effectRec, effectAdd, stdErrRec, stdErrDom, stdErrAdd, weight,row_list)


                    if robust_method == 'MERT':
                        # print(robust_method)
                         corrRecDom,corrRecAdd,corrDomAdd = cont_robust.correlations(row_list,stdErrRec, stdErrDom, stdErrAdd)
                         effect_size = cont_robust.ContRobustMert(effectDom, effectRec, corrRecDom, weight)

                # if inheritance_model != 'all':
                #     weight = 1 / var

                # get effect_sizes and variances from the model
                es_list.append(effect_size)
                var_list.append(var)
              

                n_i_list.append(N)

            bayesian_result = bayesian_meta_analysis_discrete(es_list,var_list,n_i_list)

            bayesian_result.insert(0, snp_name)# insert snp
            bayesian_result.insert(1, chrom[0]) # insert chrom
            bayesian_result.insert(2, int(pos[0]) ) # insert bp

     

            bayesian_results_list.append( bayesian_result)
    return pd.DataFrame(bayesian_results_list, columns=['SNP', 'CHR', 'BP', 'P',"E_m","E_tau_square","V_mu","V_tau_square","CI_low","CI_up",'Z'])


 