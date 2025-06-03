import numpy as np
from scipy.stats import norm


#
# vars1 = [0.16500000,0.30882353,0.07386364,0.16250000,0.08382353,0.13742690,0.09383524,0.20370370,0.40476190,0.16988417]
# vars2 = [1.04000000,0.42424242,0.09475806,0.34583333,1.00729927,0.20833333,2.08000000,0.55555555,0.36781609,0.18831169]
#
#
#
# ys1 = [1.1394343,1.4469190,1.7047481,0.4700036,0.8556661,1.4403616,0.1865860,1.5040774,1.5404450,1.6650078]
# ys2 = [3.2188758,1.2992830,0.6613985,3.2834143,4.9199809,1.3862944, 3.2188758,2.1972246, 2.2686835,-1.1451323]


def MMoM(ys1,ys2 ,vars1,vars2):
    ys1 = np.array(ys1)
    ys2 = np.array(ys2)
    vars1 = np.array(vars1)
    vars2 = np.array(vars2)

    w1 = 1/(vars1)
    w2 = 1/(vars2)
    y1_weight = sum(w1*ys1)/sum(w1)
    y2_weight = sum(w2*ys2)/sum(w2)
    # Conver the boolean list to list
    vars1_bool = (vars1 > 10**4)
    vars2_bool = (vars2 > 10**4)

    vars1_final = [i for i, x in enumerate(vars1_bool) if x]
    vars2_final = [i for i, x in enumerate(vars2_bool) if x]
    n1 = 0
    n2 = 0
    for i in range (len(vars1_final)):
        n1 += 1-1*vars1_final[i]

    for i in range (len(vars2_final)):
        n2 += 1-1*vars2_final[i]
    Q1 = sum(w1*(ys1-y1_weight)**2)
    Q2 = sum(w2*(ys2-y2_weight)**2)
    tau1_2_hat = max(0, (Q1-(n1-1))/(sum(w1)-sum(w1**2)/sum(w1)))
    tau2_2_hat = max(0, (Q2-(n2-1))/(sum(w2)-sum(w2**2)/sum(w2)))
    w1_star = 1/(vars1 + tau1_2_hat)
    w2_star = 1/(vars2 + tau2_2_hat)
    beta1_hat = sum(w1_star*ys1)/sum(w1_star)
    beta2_hat = sum(w2_star*ys2)/sum(w2_star)
    var_beta1_hat = 1/sum(w1_star)
    var_beta2_hat = 1/sum(w2_star)
    mycov_beta = sum((w1_star/sum(w1_star))*(w2_star/sum(w2_star))*(ys1 - beta1_hat)*(ys2 - beta2_hat))
    beta_hat = np.array([beta1_hat,beta2_hat])
    sigma_hat = np.matrix([[var_beta1_hat,mycov_beta],
                           [mycov_beta,var_beta2_hat]], dtype = float)
    result = {"beta.hat":beta_hat,"beta.cov":sigma_hat}
    return result


def MMoM_analysis(ys1,ys2,vars1,vars2):
    MMoM_fit = MMoM(ys1,ys2,vars1,vars2)
    # print(MMoM_fit)

    myV = np.array([1, -1])
    delta_hat = (np.dot(myV, MMoM_fit["beta.hat"]))
    delta_se = float((np.sqrt(np.dot(np.dot(myV, MMoM_fit["beta.cov"]), myV))))
    # print(delta_hat)
    # print(delta_se)

    # obtain the estimate and standard error of beta.average = (beta1+beta2)/2
    myV2 = np.array([0.5, 0.5])

    beta_average = np.dot(myV2, MMoM_fit["beta.hat"])
    beta_average_se = float(np.sqrt(np.dot(np.dot(np.transpose(myV2),MMoM_fit["beta.cov"]), myV2)))
    z_b = beta_average/beta_average_se

    p_value_beta = 1 - norm.cdf(abs(z_b))

    z_d = delta_hat / delta_se

    p_value_delta = 1 - norm.cdf(abs(z_d))
    # print(beta_average)
    # print(beta_average_se)
    return beta_average,beta_average_se,p_value_beta,delta_hat,delta_se,p_value_delta

# print(MMoM_analysis(ys1,ys2,vars1,vars2))
