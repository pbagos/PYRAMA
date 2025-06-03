#Import
import numpy as np
import pandas as pd
import cont_model,cont_robust
import math
from scipy.stats import norm, chi2

global  weight
# Open the input file
results_list = []
# path = 'demo_user/cont_data/cont_data.txt'
#
# robust_method = 'MAX'
# type_of_effect = 'RANDOM'
# inheritance_model = 'ALL'
def meta_analysis(data,inheritance_model,robust_method,type_of_effect):

 
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

                #Check the inhertance model

                if (inheritance_model == 'ALL'):
                    effect_size_dom, var_dom = cont_model.cont_model_dominant(row_list)
                    effect_size_add, var_add = cont_model.cont_model_additive(row_list)
                    effect_size_rec, var_rec = cont_model.cont_model_recessive(row_list)


                    stdErrRec = math.sqrt(var_rec)
                    stdErrDom = math.sqrt(var_dom)
                    stdErrAdd = math.sqrt(var_add)

                    if (math.sqrt(var_dom)) == 0:
                        effectDom = 0
                    else :effectDom = effect_size_dom / math.sqrt(var_dom)
                    if (math.sqrt(var_rec)) == 0:
                        effectRec = 0
                    else: effectRec = effect_size_rec / math.sqrt(var_rec)

                    if (math.sqrt(var_add)) == 0:
                        effectAdd = 0
                    else: effectAdd = effect_size_add / math.sqrt(var_add)
                    # print(f'Effect additive {effectAdd}')
                    weight = cont_robust.robustWeight(row_list, tSquared=0)



                    # Proceed with the Robust Methods
                    if robust_method == 'MAX':
                        effect_size = cont_robust.ContRobustMax(effectDom, effectRec, effectAdd, stdErrRec, stdErrDom, stdErrAdd, weight,row_list)


                    # if robust_method == 'MIN':

                    if robust_method == 'MERT':
                            corrRecDom,corrRecAdd,corrDomAdd = cont_robust.correlations(row_list,stdErrRec, stdErrDom, stdErrAdd)
                            effect_size = cont_robust.ContRobustMert(effectDom, effectRec, corrRecDom, weight)



                if (inheritance_model == 'ADDITIVE'):
                    effect_size,var = cont_model.cont_model_additive(row_list)

                if (inheritance_model == 'RECESSIVE'):
                    effect_size,var = cont_model.cont_model_recessive(row_list)

                if (inheritance_model == 'DOMINANT'):
                    effect_size,var = cont_model.cont_model_dominant(row_list)


                 # get effect_sizes and variances from the model
                es_list.append(effect_size)
                var_list.append(var)
                if (inheritance_model != 'ALL'):
                    weight = 1/var

                 

                weight_list.append(weight)
                # If Robust methods are not selected (add if statements if model = 'ALL')
                sumW += weight
                product = effect_size * weight
                sumWY += product
                sumWYSquared += effect_size * product
                sumWSquared += weight * weight
                product = 1.96 * math.sqrt(1 / weight)
                LowerLimit = effect_size - product
                UpperLimit = effect_size + product

            # print(np.array(pd.DataFrame(snp_hash['rs111'])[0]))
            # print((sumWY * sumWY / sumW))
            tSquared = 0.0

            totalVariance = sumWYSquared - (sumWY * sumWY) / sumW
            # print(f'DOF {len(XAA)}')
            degreesOfFreedom = len(XAA) - 1  # length of studies - 1
            if (totalVariance == 0):
                heterogeneity = 'NA'
            else:  heterogeneity = ((totalVariance - degreesOfFreedom) / totalVariance) * 100

            if (type_of_effect == 'RANDOM'):
                if (totalVariance > degreesOfFreedom):
                    scalingFactor = sumW - (sumWSquared / sumW)
                    tSquared = (totalVariance - degreesOfFreedom) / scalingFactor

                sumW = 0.0
                sumWY = 0.0
                # for each study the SNP exists...
                for i in range(len(XAA)):
                    if (inheritance_model == 'ALL'):

                        new_weight = cont_robust.robustWeight(row_list, tSquared)
                    else:
                        new_weight = 1 / (1 / (weight_list[i] + tSquared))
                sumW += new_weight
                sumWY += es_list[i] * new_weight
                product = 1.96 * math.sqrt(1 / new_weight)
                LowerLimit = es_list[i] - product
                UpperLimit = es_list[i] + product
            weightedMean = sumWY / sumW
            product = 1.96 * math.sqrt(1 / sumW)

            normalDib = norm.cdf(abs(weightedMean / math.sqrt(1 / sumW)))
            pValue = 2 * (1 - normalDib)
            
            # Assuming weightedMean and sumW are already defined
            z_score = abs(weightedMean / math.sqrt(1 / sumW))
            pValue = 2 * norm.sf(z_score)
            

            result = [snp_name,chrom[0],pos[0], pValue,z_score ,weightedMean - product, weightedMean + product, weightedMean]

            results_list.append(result)
            # print(result)
        # print(pd.DataFrame(results_list,columns = ['SNP','chr','p_value','Weighted_Mean - product','Weighted_Mean + product','Heterogeneity','Weighted_Mean']))
    return pd.DataFrame(results_list,columns = ['SNP','CHR','BP','P','Z','CI_lower','CI_upper','Weighted_Mean'])






