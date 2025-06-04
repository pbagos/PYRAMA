import numpy as np
from scipy.stats import norm
import pandas as pd
import numpy as np
import re
import sys
import argparse
from scipy.stats import chi2,norm

# #
# # vars1 = [0.16500000,0.30882353,0.07386364,0.16250000,0.08382353,0.13742690,0.09383524,0.20370370,0.40476190,0.16988417]
# # vars2 = [1.04000000,0.42424242,0.09475806,0.34583333,1.00729927,0.20833333,2.08000000,0.55555555,0.36781609,0.18831169]
# #
# #
# #
# # ys1 = [1.1394343,1.4469190,1.7047481,0.4700036,0.8556661,1.4403616,0.1865860,1.5040774,1.5404450,1.6650078]
# # ys2 = [3.2188758,1.2992830,0.6613985,3.2834143,4.9199809,1.3862944, 3.2188758,2.1972246, 2.2686835,-1.1451323]


# def MMoM(ys1,ys2 ,vars1,vars2):
    # ys1 = np.array(ys1)
    # ys2 = np.array(ys2)
    # vars1 = np.array(vars1)
    # vars2 = np.array(vars2)

    # w1 = 1/(vars1)
    # w2 = 1/(vars2)
    # y1_weight = sum(w1*ys1)/sum(w1)
    # y2_weight = sum(w2*ys2)/sum(w2)
    # # Conver the boolean list to list
    # vars1_bool = (vars1 > 10**4)
    # vars2_bool = (vars2 > 10**4)

    # vars1_final = [i for i, x in enumerate(vars1_bool) if x]
    # vars2_final = [i for i, x in enumerate(vars2_bool) if x]
    # n1 = 0
    # n2 = 0
    # for i in range (len(vars1_final)):
        # n1 += 1-1*vars1_final[i]

    # for i in range (len(vars2_final)):
        # n2 += 1-1*vars2_final[i]
    # Q1 = sum(w1*(ys1-y1_weight)**2)
    # Q2 = sum(w2*(ys2-y2_weight)**2)
    # tau1_2_hat = max(0, (Q1-(n1-1))/(sum(w1)-sum(w1**2)/sum(w1)))
    # tau2_2_hat = max(0, (Q2-(n2-1))/(sum(w2)-sum(w2**2)/sum(w2)))
    # w1_star = 1/(vars1 + tau1_2_hat)
    # w2_star = 1/(vars2 + tau2_2_hat)
    # beta1_hat = sum(w1_star*ys1)/sum(w1_star)
    # beta2_hat = sum(w2_star*ys2)/sum(w2_star)
    # var_beta1_hat = 1/sum(w1_star)
    # var_beta2_hat = 1/sum(w2_star)
    # mycov_beta = sum((w1_star/sum(w1_star))*(w2_star/sum(w2_star))*(ys1 - beta1_hat)*(ys2 - beta2_hat))
    # beta_hat = np.array([beta1_hat,beta2_hat])
    # sigma_hat = np.matrix([[var_beta1_hat,mycov_beta],
                           # [mycov_beta,var_beta2_hat]], dtype = float)
    # result = {"beta.hat":beta_hat,"beta.cov":sigma_hat}
    # return result


# def MMoM_analysis(ys1,ys2,vars1,vars2):
    # MMoM_fit = MMoM(ys1,ys2,vars1,vars2)
    # # print(MMoM_fit)

    # myV = np.array([1, -1])
    # delta_hat = (np.dot(myV, MMoM_fit["beta.hat"]))
    # delta_se = float((np.sqrt(np.dot(np.dot(myV, MMoM_fit["beta.cov"]), myV))))
    # # print(delta_hat)
    # # print(delta_se)

    # # obtain the estimate and standard error of beta.average = (beta1+beta2)/2
    # myV2 = np.array([0.5, 0.5])

    # beta_average = np.dot(myV2, MMoM_fit["beta.hat"])
    # beta_average_se = float(np.sqrt(np.dot(np.dot(np.transpose(myV2),MMoM_fit["beta.cov"]), myV2)))
    # z_b = beta_average/beta_average_se

    # p_value_beta = 1 - norm.cdf(abs(z_b))

    # z_d = delta_hat / delta_se

    # p_value_delta = 1 - norm.cdf(abs(z_d))
    # # print(beta_average)
    # # print(beta_average_se)
    # return beta_average,beta_average_se,p_value_beta,delta_hat,delta_se,p_value_delta

# # print(MMoM_analysis(ys1,ys2,vars1,vars2))
#!/usr/bin/env python3


def mmom_multi(ys, vars_):
 
    # Extract effect sizes
    ys1 = ys[:, 0]
    ys2 = ys[:, 1]

    # Extract variances
    vars1 = vars_[:, 0]
    vars2 = vars_[:, 1]

    # Inverse-variance weights
    w1 = 1.0 / vars1
    w2 = 1.0 / vars2

    # Weighted means
    y1_weight = np.sum(w1 * ys1) / np.sum(w1)
    y2_weight = np.sum(w2 * ys2) / np.sum(w2)

    # Count number of studies not imputed (i.e., not missing)
    n1 = np.sum(vars1 < 1e4)
    n2 = np.sum(vars2 < 1e4)

    # Q statistics
    Q1 = np.sum(w1 * (ys1 - y1_weight) ** 2)
    Q2 = np.sum(w2 * (ys2 - y2_weight) ** 2)

    
        
    # Between-study variances (tau^2 estimates)
    tau1_2_hat = max(0, (Q1 - (n1 - 1)) / (np.sum(w1) - np.sum(w1**2) / np.sum(w1)))
    tau2_2_hat = max(0, (Q2 - (n2 - 1)) / (np.sum(w2) - np.sum(w2**2) / np.sum(w2)))

    # Adjusted weights
    w1_star = 1.0 / (vars1 + tau1_2_hat)
    w2_star = 1.0 / (vars2 + tau2_2_hat)

    # Final beta estimates
    beta1_hat = np.sum(w1_star * ys1) / np.sum(w1_star)
    beta2_hat = np.sum(w2_star * ys2) / np.sum(w2_star)

    # Variance of beta estimates
    var_beta1_hat = 1.0 / np.sum(w1_star)
    var_beta2_hat = 1.0 / np.sum(w2_star)

    # Covariance between beta1 and beta2
    w1_norm = w1_star / np.sum(w1_star)
    w2_norm = w2_star / np.sum(w2_star)
    mycov_beta = np.sum(w1_norm * w2_norm * (ys1 - beta1_hat) * (ys2 - beta2_hat))

    # Combine results
    beta_hat = np.array([beta1_hat, beta2_hat])
    sigma_hat = np.array([
        [var_beta1_hat, mycov_beta],
        [mycov_beta, var_beta2_hat]
    ])

    return {
        'beta_hat': beta_hat,
        'beta_cov': sigma_hat
    }


def beta_SE_meta(data):
  
    df =  data.dropna()
    
    # Order by the 'SNP' column.
    df = df.sort_values(by='SNP')
    
    # Identify BETA and SE columns using regex (e.g., BETA1, BETA2, ... and SE1, SE2, ...).
    beta_cols = sorted([col for col in df.columns if re.match(r'BETA\d+', col)],
                       key=lambda x: int(re.findall(r'\d+', x)[0]))
    se_cols = sorted([col for col in df.columns if re.match(r'SE\d+', col)],
                     key=lambda x: int(re.findall(r'\d+', x)[0]))
    
    # Check that each BETA column has a matching SE column.
    if len(beta_cols) != len(se_cols):
        print("Error: The number of BETA columns does not match the number of SE columns.")
        sys.exit(1)
    
    results = []
    # Group the dataframe by 'SNP'.
    grouped = df.groupby('SNP')
    
    for var, group in grouped:
        row = {'SNP': var}
        
        # Build matrices for effect estimates (ys) and standard errors.
        # Each column corresponds to one method.
        ys_mat = np.array(np.column_stack([group[col].values for col in beta_cols]),dtype = np.float64 )
        se_mat = np.array(np.column_stack([group[col].values for col in se_cols]),dtype = np.float64  )
        # Convert standard errors to variances.
        vars_mat = np.array(se_mat,dtype = np.float64) ** 2
        
        # Run the multivariate meta-analysis (MMoM_multi equivalent).
        res = mmom_multi(ys_mat, vars_mat)
        beta_hat = res['beta_hat']
        beta_cov = res['beta_cov']
        #print(beta_hat)
       # print(beta_cov)
        # Save results for each method: meta_BETA and meta_SE (sqrt of variance).
        for k, (bh, var_bh) in enumerate(zip(beta_hat, np.diag(beta_cov)), start=1):
            row[f'meta_BETA{k}'] = bh
            row[f'meta_SE{k}'] = np.sqrt(var_bh)
        
        # Compute the Wald statistic and associated p-value.
        try:
            inv_cov = np.linalg.pinv(beta_cov)
            # Wald statistic: beta_hat' * inv_cov * beta_hat.
            wald = beta_hat.T @ inv_cov @ beta_hat
            p_value = chi2.sf(wald, df=beta_cov.shape[0])
        except np.linalg.LinAlgError:
            wald = np.nan
            p_value = np.nan
        
        row["Wald"] = wald
        row["P"] = p_value
        
        results.append(row)
    
    # Create a DataFrame for the results and save it to a tab-separated file.
    out_df = pd.DataFrame(results)
    return out_df
 


def main(input_file, output_file):
    import time 
    start= time.time()
    # Read the input tab-separated file.
    df = pd.read_csv(input_file, sep="\t")
    
    # Order by the 'SNP' column.
    df = df.sort_values(by='SNP')
    
    # Identify BETA and SE columns using regex (e.g., BETA1, BETA2, ... and SE1, SE2, ...).
    beta_cols = sorted([col for col in df.columns if re.match(r'BETA\d+', col)],
                       key=lambda x: int(re.findall(r'\d+', x)[0]))
    se_cols = sorted([col for col in df.columns if re.match(r'SE\d+', col)],
                     key=lambda x: int(re.findall(r'\d+', x)[0]))
    
    # Check that each BETA column has a matching SE column.
    if len(beta_cols) != len(se_cols):
        print("Error: The number of BETA columns does not match the number of SE columns.")
        sys.exit(1)
    
    results = []
    # Group the dataframe by 'SNP'.
    grouped = df.groupby('SNP')
    
    for var, group in grouped:
        row = {'SNP': var}
        
        # Build matrices for effect estimates (ys) and standard errors.
        # Each column corresponds to one method.
        ys_mat = np.column_stack([group[col].values for col in beta_cols])
        se_mat = np.column_stack([group[col].values for col in se_cols])
        # Convert standard errors to variances.
        vars_mat = se_mat ** 2
        
        # Run the multivariate meta-analysis (MMoM_multi equivalent).
        res = mmom_multi(ys_mat, vars_mat)
        beta_hat = res['beta_hat']
        beta_cov = res['beta_cov']
        #print(beta_hat)
       # print(beta_cov)
        # Save results for each method: meta_BETA and meta_SE (sqrt of variance).
        for k, (bh, var_bh) in enumerate(zip(beta_hat, np.diag(beta_cov)), start=1):
            row[f'meta_BETA{k}'] = bh
            row[f'meta_SE{k}'] = np.sqrt(var_bh)
        
        # Compute the Wald statistic and associated p-value.
        try:
            inv_cov = np.linalg.pinv(beta_cov)
            # Wald statistic: beta_hat' * inv_cov * beta_hat.
            wald = beta_hat.T @ inv_cov @ beta_hat
            p_value = chi2.sf(wald, df=beta_cov.shape[0])
        except np.linalg.LinAlgError:
            wald = np.nan
            p_value = np.nan
        
        row["Wald"] = wald
        row["P"] = p_value
        
        results.append(row)
    
    # Create a DataFrame for the results and save it to a tab-separated file.
    out_df = pd.DataFrame(results)
    end = time.time()
    print(f"Executed in : {end-start} seconds.")

    out_df.to_csv(output_file, index=False, sep="\t")
    print(f"Meta-analysis results have been saved to: {output_file}")


if __name__ == '__main__':
  
    parser = argparse.ArgumentParser(
        description='Multivariate random-effects meta-analysis for beta and SE pairs in a tab-separated file')
    parser.add_argument('--input_file',
                        help='Input tab-separated file with columns: SNP, BETA1, SE1, BETA2, SE2, ...')
    parser.add_argument('--output_file', help='Output tab-separated file to save the meta-analysis results')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
 
    