import numpy as np
import math
from scipy.integrate import quad
from scipy.stats import norm, chi2
from scipy.stats import chi2_contingency
from scipy.stats import cauchy
from decimal import Decimal, getcontext

# Set the precision 
getcontext().prec = 1000

global new_weight


def discreteCorrPrototypeFunction(a, b, c):
    numerator = a * (b + 2 * c)

    denominator = math.sqrt(a * (1 - a)) * math.sqrt((b + 2 * a) * c + (b + 2 * c) * a)

    if denominator == 0:
        return 0
    return numerator / denominator


def setProb1(prob1):
    if (prob1() >= 0.0 and prob1 <= 1.0):
        return prob1
    else:
        print("Error: Prob1 out of bounds [0-1 ]")
        exit()


def setProb2(prob2):
    if (prob2() >= 0.0) & (prob2 <= 1.0):
        return prob2
    else:
        print("Error: Prob1 out of bounds [0-1 ]")
        exit()


def setDirection(direction):
    if (direction == 1) | (direction == -1):
        return direction
    else:
        print("Error: direction must be either 1 or -1")
        exit()


def findMax(effectDom, effectRec, effectAdd):
    additive = abs(effectAdd)
    recessive = abs(effectRec)
    max = abs(effectDom)

    if (additive > max):
        max = additive
    if (recessive > max):
        max = recessive

    return max


def gaussLegendreQuad(max0, w0, w1, corr_rec_dom):
    b1 = max0 * (1 - w1) / w0

    def f1(x):
        return norm.cdf((max0 - corr_rec_dom * x) / math.sqrt(1 - corr_rec_dom * corr_rec_dom)) * norm.pdf(x)

    def f2(x):
        return norm.cdf(
            ((max0 - w0 * x) / w1 - corr_rec_dom * x) / math.sqrt(1 - corr_rec_dom * corr_rec_dom)) * norm.pdf(x)

    def f3(x):
        return norm.cdf((-max0 - corr_rec_dom * x) / math.sqrt(1 - corr_rec_dom * corr_rec_dom)) * norm.pdf(x)

    I1, r1 = quad(f1, 0, max0 * (1 - w1) / w0, limit=100)
    I2, r2 = quad(f2, max0 * (1 - w1) / w0, max0, limit=100)
    if I1 >= 0.49999999999999999:
        print(f'I1: {I1}, I2: {I2} ')
        return abs(norm.ppf((2 * I2) / 2))

    I3, r3 = quad(f3, 0, max0, limit=100)

    calculation = 2 * Decimal(I1) + 2 * Decimal(I2) - 2 * Decimal(I3)
    # Define p1 as a very small number
    calculation = np.float64(calculation)
    # #calculation = np.exp(np.log(2*I1+2*I2)) -2*I3
    # #calculation  = np.exp(np.log(2*I1) + np.log(1 + np.exp(np.log(2*I2) - np.log(2*I1)))) -2*I3
    # #calculation  = np.exp(math.log(2 * I1) + math.log(2 * I2) - math.log(2 * I3))

    # # Compute y'
    # y_prime = np.log(I1) + np.log(1 + np.exp(np.log(I2) - np.log(I3)) - np.exp(np.log(I3) - np.log(I1))) + np.exp(np.log(-I3) - np.log(I1))
    # 
    # # Compute poy
    # calculation = 1 - 2 * np.exp(y_prime)
    if calculation > 1 :
        p1 = norm.sf(max0) * 2
        p1 = Decimal(p1)

        p = 1 - ((1 - p1) ** 3)

        p = np.float64(p)
        zValue = abs(norm.ppf(p / 2) )

        return zValue
    else :
        return abs(norm.ppf( (1 - calculation) / 2))


def DiscteteRobustMert(effectDom, effectRec, corrRecDom, weight):
    return ((effectDom + effectRec) / math.sqrt(2 * (1 + corrRecDom))) / math.sqrt(weight)


def DiscreteRobustMax(effectDom, effectRec, effectAdd, corrDomAdd, corrRecAdd, corrRecDom, robustWeight,
                      approximate_max):
    sign = np.sign(effectAdd)
    max_ = max(abs(effectDom), abs(effectRec), abs(effectAdd))

    w0 = (corrRecAdd - corrRecDom * corrDomAdd) / (1 - corrRecDom * corrRecDom)
    w1 = (corrDomAdd - corrRecDom * corrRecAdd) / (1 - corrRecDom * corrRecDom)

    if approximate_max == 'YES':
        # #Simple Case
        # zValue = abs(norm.ppf((norm.sf(abs(max_) * 2 / 2.05))))

        #
        #Cauchy Method
        # list_p  = [norm.sf(abs(effectDom ) )  * 2,norm.sf(abs(effectRec)) * 2,norm.sf(abs(effectAdd * 2)) * 2]
        # T = ((1/3) *  np.tan( (0.5 - list_p[0])*np.pi ) )  +    ((1/3) *  np.tan( (0.5 - list_p[1])*np.pi ) )       +((1/3) *  np.tan( (0.5 - list_p[2])*np.pi ) )
        # pValue = cauchy.sf(abs( T )) * 2
        # zValue = norm .ppf(pValue/2)

        #MinP Method

        p1 = norm.sf(max_)*2
        # print(p1)
        # Define p1 as a very small number
        p1 = Decimal(p1)

        # Perform the calculation using Decimal
        p = 1 - ((1 - p1) ** 3)

        # Output the result
        #print("p =", p)

        #p = 1 - ((1 - p1) ** 3)
        # if p == 0 :
        #     zValue = abs(norm.ppf((norm.sf(abs(max_) * 2 / 2.05))))
        # else:
        p = np.float64(p)
        zValue = norm.ppf(p/2)





    else:
       zValue = gaussLegendreQuad(max_, w0, w1, corrRecDom)
    EffectSize = ((zValue * sign) / math.sqrt(robustWeight))
    # if abs(EffectSize) == np.inf :
    #     print (zValue,max_,w0,w1,corrRecDom)
    # print('-------------------------')
    # # print(zValue,EffectSize)
    # print('-------------------------')

    return EffectSize


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
    # print(f"Observed genotypes {observed}")

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


def DiscreteRobustMin(weight, effectAdd, pChiSquared):
    sign = math.copysign(1, effectAdd)
    # print('Sign' + str(sign))
    pAdd = 2 * norm.cdf(-abs(effectAdd))
    # pAdd = 2 * norm.cdf(abs(effectAdd))

    # print(-abs(effectAdd))
    # print(pAdd)
    min_ = min(pAdd, pChiSquared)
    # print('pchi2 '+str(pChiSquared))
    # print(min_)
    if pChiSquared < pAdd or pChiSquared == pAdd:
        min_ = pChiSquared
    else:
        min_ = pAdd
    if min_ == 0:
        min_ = 1e-285

    tAdd = norm.ppf(min_ / 2) * norm.ppf(min_ / 2)

    tMethod = -2 * math.log(min_)

    # else:
    #     tMethod = 0

    # print(f"Tadd {tAdd} and  tMethod {tMethod}")
    if tAdd < tMethod:
        pMin = 0.5 * math.exp(-tAdd / 2) + 0.5 * (min_) - (1 / (2 * math.pi)) * gauss_legendre_quad(tMethod, tAdd)
        # pMin = 0.5 * math.exp(-tAdd / 2) - 0.5 * math.exp(-min_/2) + (1 / (2 * math.pi)) * gauss_legendre_quad(min_, tAdd)

        # pMin = 2.2 * min_
    else:
        pMin = min_


    zValue = abs(norm.ppf(pMin / 2))
    # MINP_ method
    #zValue = abs(norm.ppf((1-(1-min_)**2) / 2))
    effect_size = zValue * sign / math.sqrt(weight)
    return effect_size


def gauss_legendre_quad(tMethod, tAdd):
    def f1(x):
        return math.exp(-x / 2) * np.arcsin(((2 * tAdd) / x) - 1)

    I1, r1 = quad(f1, tAdd, tMethod, limit=100, points=1000)

    # print(I1)
    return I1


def robustWeight(R, S, tSquared):
    inverse = 1 / (1 / R + 1 / S)
    weight = 1 / (1 / inverse + tSquared)
    return weight


def model_dominant(row_list, effect_size_type):
    effect_size = 0
    var = 0
    aa1 = float(row_list[0])
    ab1 = float(row_list[1])
    bb1 = float(row_list[2])
    aa0 = float(row_list[3])
    ab0 = float(row_list[4])
    bb0 = float(row_list[5])
    R = aa1 + ab1 + bb1
    S = aa0 + ab0 + bb0
    N = R + S
    n2 = bb0 + bb1
    n1 = ab1 + ab0
    if (effect_size_type == 'OR'):
        effect_size = math.log(((aa1 + ab1) / (aa0 + ab0)) / (bb1 / bb0))
        var = 1 / (aa1 + ab1) + 1 / (aa0 + ab0) + 1 / bb1 + 1 / bb0

    # if (effect_size_type == 'RR'):
    #     effect_size = ((aa1 + ab1) * (bb1 + bb0)) / (bb1 * (aa1 + ab1 + aa0 + ab0))
    #
    #     var = 1 / (aa1 + ab1) - 1 / (aa1 + ab1 + aa0 + ab0) + 1 / bb1 - 1 / bb1 / bb0
    #
    # if (effect_size_type == 'RISK_DIFF'):
    #     effect_size = (aa1 + ab1) / (aa1 + ab1 + aa0 + ab0) - bb1 / (bb1 + bb0)
    #
    #     var = ((aa1 + ab1) * (aa0 + ab0)) / (aa1 + ab1 + aa0 + ab0) ** 3 - bb1 / (bb1 + bb0)

    if (effect_size_type == 'CATT'):
        effect_size = 1 / N * ((S * (ab1 + bb1)) - (R * (ab0 + bb0)))

        var = (R * S / N) * (((n1) + n2) / N - (((n1 + n2) / N) * ((n1 + n2) / N)))

    return effect_size, var


def model_additive(row_list, effect_size_type):  # allelic model
    effect_size = 0
    var = 0
    aa1 = float(row_list[0])
    ab1 = float(row_list[1])
    bb1 = float(row_list[2])
    aa0 = float(row_list[3])
    ab0 = float(row_list[4])
    bb0 = float(row_list[5])
    R = aa1 + ab1 + bb1
    S = aa0 + ab0 + bb0
    N = R + S
    n2 = bb0 + bb1
    n1 = ab1 + ab0

    if (effect_size_type == 'OR'):
        effect_size = math.log(((2 * bb1 + ab1) / (2 * bb0 + ab0)) / ((2 * aa1 + ab1) / (2 * aa0 + ab0)))
        var = 1 / (2 * bb1 + ab1) + 1 / (2 * bb0 + ab0) + 1 / (2 * aa1 + ab1) + 1 / (2 * aa0 + ab0)

    # if (effect_size_type == 'RR'):
    #     effect_size = math.log(((2 * bb1 + ab1) * ((2 * aa1 + ab1) + (2 * aa0 + ab0))) / (
    #             (2 * aa1 + ab1) * (2 * bb1 + ab1 + 2 * bb0 + ab0)))
    #
    #     var = 1 / (2 * bb1 + ab1) - 1 / ((2 * bb1 + ab1) + (2 * bb0 + ab0)) + 1 / (2 * aa1 + ab1) - 1 / (
    #             2 * aa1 + ab1 + 2 * aa0 + ab0)

    # if (effect_size_type == 'RISK_DIFF'):
    #     effect_size = (2 * bb1 + ab1) / (2 * bb1 + ab1 + 2 * bb0 + ab0) - (2 * aa1 + ab1) / (
    #             2 * aa1 + ab1 + 2 * aa0 + ab0)
    #
    #     var = ((2 * bb1 + ab1) * (2 * bb0 + ab0)) / (2 * bb1 + ab1 + 2 * bb0 + ab0) ** 3 + (
    #             (2 * aa1 + ab1) * (2 * aa0 + ab0)) / (2 * aa1 + ab1 + 2 * aa0 + ab0) ** 3

    if (effect_size_type == 'CATT'):
        effect_size = 1 / N * ((S * (0.5 * ab1 + bb1)) - (R * (0.5 * ab0 + bb0)))

        var = (R * S / N) * (((0.5 * 0.5 * n1) + n2) / N - (((0.5 * n1 + n2) / N) * ((0.5 * n1 + n2) / N)))

    return effect_size, var


def model_recessive(row_list, effect_size_type):
    effect_size = 0
    var = 0
    aa1 = float(row_list[0])
    ab1 = float(row_list[1])
    bb1 = float(row_list[2])
    aa0 = float(row_list[3])
    ab0 = float(row_list[4])
    bb0 = float(row_list[5])
    R = aa1 + ab1 + bb1
    S = aa0 + ab0 + bb0
    N = R + S
    n2 = bb0 + bb1
    n1 = ab1 + ab0

    # print(aa1,ab1,bb1,aa0,ab0,bb0)
    if (effect_size_type == 'OR'):
        effect_size = math.log(((bb1 + ab1) /
                                (bb0 + ab0)) /
                               (aa1 / aa0))

        var = 1 / (bb1 + ab1) + 1 / (bb0 + ab0) + 1 / aa1 + 1 / aa0
    #
    # if (effect_size_type == 'RR'):
    #     effect_size = ((aa1 + ab1) * (bb1 + bb0)) / (bb1 * (aa1 + ab1 + ab0 + ab0))
    #     var = 1 / (aa1 + ab1) - 1 / (aa1 + ab1 + aa0 + ab0) + 1 / bb1 - 1 / (bb1 + bb0)
    #
    # if (effect_size_type == 'RISK_DIFF'):
    #     effect_size = (aa1 + ab1) / (aa1 + ab1 + aa0 + ab0) - bb1 / (bb1 + bb0)
    #
    #     var = ((aa1 + ab1) * (aa0 + ab0)) / (aa1 + ab1 + aa0 + ab0) ** 3 + (bb1 * bb0) / (bb1 + bb0) ** 3

    if (effect_size_type == 'CATT'):
        effect_size = 1 / N * (S * bb1 - R * bb0)

        var = (R * S / N) * ((n2 / N) - ((n2 / N) * (n2 / N)))
    return effect_size, var
