import numpy as np
import math
import sys
from scipy.integrate import quad
from scipy.stats import norm, chi2
from scipy.stats import chi2_contingency


def robustWeight(row_list,tSquared):
    xaa = row_list[0]
    sdaa = row_list[1]
    naa = row_list[2]
    xab = row_list[3]
    sdab = row_list[4]
    nab = row_list[5]
    xbb = row_list[6]
    sdbb = row_list[7]
    nbb = row_list[8]

    totalPopulation = naa +nbb + nab
    p = 2 * naa + nab
    q = 2 * nbb + nab
    phenotypeVariance = ((naa - 1) * (sdaa * sdaa) + (nbb - 1) *(sdbb * sdbb) + (nab - 1) * (sdab * sdab))/ (nbb + naa + nab - 3)
    
    
    phenotypeVariance = (
        ((naa - 1) * (sdaa ** 2)) +
        ((nbb - 1) * (sdbb ** 2)) +
        ((nab - 1) * (sdab ** 2))
    ) / (naa + nbb + nab - 3)
    
    # print(f'Total Population {totalPopulation}')
    # print(f'p {p}')
    # print(f'q {q}')
    # print(f' phenotypeVariance {phenotypeVariance}')
    Ib = (2 * totalPopulation * p * q) / phenotypeVariance
    # print(Ib)
    weight = 1 / (1 / Ib + tSquared)
    # print(weight)
    return weight

def gaussLegendreQuad(max0, w0, w1, corr_rec_dom):
    b1 = max0 * (1 - w1) / w0

    def f1(x):
        return norm.cdf((max0 - corr_rec_dom * x) / math.sqrt(1 - corr_rec_dom * corr_rec_dom)) * norm.pdf(x)

    def f2(x):
        return norm.cdf(
            ((max0 - w0 * x) / w1 - corr_rec_dom * x) / math.sqrt(1 - corr_rec_dom * corr_rec_dom)) * norm.pdf(x)

    def f3(x):
        return norm.cdf((-max0 - corr_rec_dom * x) / math.sqrt(1 - corr_rec_dom * corr_rec_dom)) * norm.pdf(x)

    I1, r1 = quad(f1, 0, max0 * (1 - w1) / w0)
    I2, r2 = quad(f2, max0 * (1 - w1) / w0, max0)
    I3, r3 = quad(f3, 0, max0)

    calculation = 2 * I1 + 2 * I2 - 2 * I3

    # print("Calculation: " + calculation)
    return abs(norm.ppf((1 - calculation) / 2))


def correlations(row_list, stdErrRec, stdErrDom, stdErrAdd):
    xaa = row_list[0]
    sdaa = row_list[1]
    naa = row_list[2]
    xab = row_list[3]
    sdab = row_list[4]
    nab = row_list[5]
    xbb = row_list[6]
    sdbb = row_list[7]
    nbb = row_list[8]
    
    covDomAdd = ((2 * nbb * sdbb * sdbb) + (nab * sdab * sdab)) \
            / ((nbb + nab) * (2 * nbb + nab)) \
            + (2 * sdaa * sdaa * (nbb + nab) - (nab * sdab * sdab)) \
            / ((nbb + nab) * (2 * naa  + nab))
    covRecDom = (sdbb * sdbb * (nab + naa) + (sdaa * sdaa) * (nab + nbb)
            - (sdab * sdab * nab)) \
            / ((nab + naa) * (nab + nbb))
    covRecAdd = ((2 * naa * sdaa * sdaa) + (nab * sdab * sdab)) \
           / ((naa + nab) * (2 * naa + nab)) \
            + (2 * sdbb * sdbb * (naa + nab) - (nab * sdab * sdab)) \
           / ((naa + nab) * (2 * nbb  + nab))
    corrRecDom = covRecDom / (stdErrDom * stdErrRec)
    corrRecAdd = covRecAdd / (stdErrAdd * stdErrRec)
    corrDomAdd = covDomAdd / (stdErrDom * stdErrAdd)

    return (corrRecDom,corrRecAdd,corrDomAdd)



def ContRobustMax(effectDom, effectRec, effectAdd, stdErrRec, stdErrDom, stdErrAdd, robustWeight,row_list):
    sign = np.sign(effectAdd)
    max_ = max(abs(effectDom), abs(effectRec), abs(effectAdd))

    corrRecDom,corrRecAdd,corrDomAdd = correlations(row_list,stdErrRec, stdErrDom, stdErrAdd)
    w0 = (corrRecAdd - corrRecDom * corrDomAdd) / (1 - corrRecDom * corrRecDom)
    w1 = (corrDomAdd - corrRecDom * corrRecAdd) / (1 - corrRecDom * corrRecDom)
    zValue = gaussLegendreQuad(max_, w0, w1, corrRecDom)
    EffectSize = ((zValue * sign) / math.sqrt(robustWeight))

    return EffectSize

def ContRobustMert(effectDom, effectRec, corrRecDom, weight):
    return ((effectDom + effectRec) / math.sqrt(2 * (1 + corrRecDom))) / math.sqrt(weight)


def pValueChiSquared(row_list):
    aa1 = row_list[0]
    ab1 = row_list[1]
    bb1 = row_list[2]
    aa0 = row_list[3]
    ab0 = row_list[4]
    bb0 = row_list[5]
    R = aa1 + ab1 + bb1  # cases
    S = aa0 + ab0 + bb0  # controls
    N = R + S

    # Create 2x3 contingency table
    observed = np.array([[aa1, ab1, bb1],
                         [aa0, ab0, bb0]])
    print(f"Observed genotypes {observed}")

    # Calculate chi-squared statistic, p-value, degrees of freedom, and expected frequencies
    chi2_stat, p_val, dof, expected = chi2_contingency(observed)




    # observed = [aa0, ab0, bb0,
    #             aa1, ab1, bb1]
    # genotypes = [aa1 + aa0,
    #              ab1 + ab0,
    #              bb1 + bb0]
    # expected = [0] * 6
    # genLength = len(genotypes)
    # cases = R
    # controls = S
    # total = N
    # xSquared = 0
    # for i in range(len(observed)):
    #     if i < 3:
    #         expected[i] = (genotypes[i] * controls) / total
    #     else:
    #         expected[i] = (genotypes[i % genLength] * cases) / total
    #     if expected[i] ==0 :
    #         xSquared += 0
    #     else: xSquared += ((observed[i] - expected[i]) * (observed[i] - expected[i])) / expected[i]
    return p_val


def ContRobustMin(weight, effectAdd, pChiSquared):
    sign = math.copysign(1, effectAdd)
    print('Sign' + str(sign))
    pAdd = 2 * norm.cdf(-abs(effectAdd))
    # pAdd = 2 * norm.cdf(abs(effectAdd))

    print(-abs(effectAdd))
    print(pAdd)
    min_ = min(pAdd, pChiSquared)
    print('pchi2 '+str(pChiSquared))
    print(min_)
    if min_ == 0:
        min_ = 1e-100

    tAdd =  norm.ppf(min_ / 2) * norm.ppf(min_ / 2)

    tMethod = -2 * math.log(min_)

    # else:
    #     tMethod = 0

    print(f"Tadd {tAdd} and  tMethod {tMethod}")
    if tAdd < tMethod:
        pMin = 0.5 * math.exp(-tAdd / 2) + 0.5 * (min_) + (1 / (2 * math.pi)) * gauss_legendre_quad(tMethod, tAdd)
    else:
        pMin = min_
    zValue = abs(norm.ppf(pMin / 2))
    effect_size = zValue * sign / math.sqrt(weight)
    return effect_size


def gauss_legendre_quad(tMethod, tAdd):
    def f1(x):
        return math.exp(-x / 2) * np.arcsin(((2 * tAdd) / x) - 1)

    I1, r1 = quad(f1, tAdd, tMethod)


    print(I1)
    return I1
