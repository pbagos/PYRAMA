import os
import pandas as pd
import numpy as np
import math
import robust, discrete_model
import sys
from scipy.integrate import quad
from scipy.stats import norm, chi2
from bayes_factor import bayes_factor

global new_weight, results_list, bf

results_list = []


def meta_analysis(data,inheritance_model, effect_size_type, robust_method, type_of_effect, approximate_max):
    # study_list = []
    # for file in file_list:
    #     study_list.append(pd.read_csv(path + file, sep='\t', encoding='latin1'))
    # merge... files
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

            bf = bayes_factor(abs(effect_size), np.sqrt(var))
            # print(f'BF : {bf}')
            # print("CASES AND CONTROLS "+ str(R)+"  "+ str(S))
            # print('WEIGHT:' + str(weight))
            # print('E_SIZE:' + str(effect_size))
            if (inheritance_model == 'ALL'):
                weight = robust.robustWeight(R, S, tSquared=0)
            else:
                weight = 1 / var

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

        totalVariance = sumWYSquared - (sumWY * sumWY) / abs(sumW)
        # print(f'DOF {len(AA1)}')
        if len(AA1) > 1:
            degreesOfFreedom = len(AA1) - 1  # length of studies - 1
        else:
            degreesOfFreedom = 1

        if (totalVariance == 0):
            heterogeneity = 'NA'
        else:
            heterogeneity = ((totalVariance - degreesOfFreedom) / totalVariance) * 100

        if (type_of_effect == 'RANDOM'):
            if (totalVariance > degreesOfFreedom):

                scalingFactor = sumW - (sumWSquared / sumW)
                if scalingFactor != 0:
                    tSquared = (totalVariance - degreesOfFreedom) / scalingFactor
                else:
                    # Handle the case when scalingFactor is zero (division by zero)
                    tSquared = 0  # You can assign an appropriate default value or handle it differently

            sumW = 0.0
            sumWY = 0.0
            # for each study the SNP exists...
            for i in range(len(AA1)):
                if (inheritance_model == 'ALL'):

                    new_weight = robust.robustWeight(R, S, tSquared)
                else:
                    new_weight = 1 / (1 / (weight_list[i] + tSquared))
            sumW += new_weight
            sumWY += es_list[i] * new_weight
            product = 1.96 * math.sqrt(1 / abs(new_weight))
            LowerLimit = es_list[i] - product
            UpperLimit = es_list[i] + product
        weightedMean = sumWY / sumW
        product = 1.96 * math.sqrt(1 / abs(sumW))

        # normalDib = norm.cdf(abs(weightedMean / math.sqrt(1 / sumW)))
        # pValue = 2 * (1 - normalDib)
        z = abs(weightedMean / math.sqrt(1 / abs(sumW)))
        pValue = norm.sf(abs(z)) * 2
        # normalDib = norm.cdf(abs(weightedMean / math.sqrt(1 / sumW)))
        # pValue = 2 * (1 - normalDib)
        result = [snp_name, chromosome, pos[0], pValue, weightedMean - product, weightedMean + product,
                  weightedMean, bf]

        results_list.append(result)
        # print(result)
    # print(pd.DataFrame(results_list,columns = ['SNP','CHR','P','Weighted_Mean - product','Weighted_Mean + product','Heterogeneity','Weighted_Mean']))
    return pd.DataFrame(results_list, columns=['SNP', 'CHR', 'BP', 'P', 'CI_lower', 'CI_upper' ,
                                               'Weighted_Mean', 'Bayes_Factor'])
