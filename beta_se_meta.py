import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats.distributions import chi2

def altmeta(y1, s2):
    # Convert variances
    variances = np.array(s2)**2
    n = len(y1)
    w = [(1 / x) for x in s2 if x != 0]  # Prevent division by zero

    if len(w) == 0 or np.isclose(sum(w), 0):
        # If all variances are zero or sum of weights is zero
        return (None, 0, 0, None, None, None, None, 0)

    # Weighted mean calculation
    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)

    # Q-statistic calculation
    Q = sum(a * ((x - mu_bar) ** 2) for a, x in zip(w, y1))
    p_Q = chi2.sf(Q, n - 1)

    # I-squared and tau-squared calculations
    if Q == 0:
        I2 = 0
        tau2_DL = 0
    else:
        I2 = max(0, (Q - (n - 1)) / Q)
        denom = sum(w) - (sum([x ** 2 for x in w]) / sum(w))
        tau2_DL = max(0, (Q - n + 1) / denom) if denom != 0 else 0

    # Adjusted weights with tau2_DL
    w = [(1 / (v + tau2_DL)) for v in variances if v + tau2_DL != 0]
    if len(w) == 0 or sum(w) == 0:
        return (None, I2, tau2_DL, p_Q, None, None, None, Q)

    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)
    se = np.sqrt(1 / sum(w))
    z = mu_bar / se if se != 0 else 0
    p = norm.sf(abs(z)) * 2

    return (p, I2, tau2_DL, p_Q, se, z, mu_bar, Q)

def beta_se(file, biv_beta_input):
    beta_hash, se_hash, beta_hash2, se_hash2 = {}, {}, {}, {}
    chrom_hash, pos_hash = {}, {}

    if biv_beta_input == 'YES':
        for _, row in file.iterrows():
            snp = row['SNP']
            beta, se, beta2, se2 = row[3], row[4], row[5], row[6]
            chrom, pos = row[1], row[2]

            # Populate dictionaries
            beta_hash.setdefault(snp, []).append(beta)
            se_hash.setdefault(snp, []).append(se)
            beta_hash2.setdefault(snp, []).append(beta2)
            se_hash2.setdefault(snp, []).append(se2)
            chrom_hash.setdefault(snp, []).append(chrom)
            pos_hash.setdefault(snp, []).append(pos)

        result_df = pd.DataFrame(columns=['SNP', 'CHR', 'BP', 'Q1', 'I2_1', 'tau2_DL_1', 'p_Q1', 'se1', 'z1', 
                                          'Weighted_Mean1', 'p_value1', 'Q2', 'I2_2', 'tau2_DL_2', 'p_Q2', 
                                          'se2', 'z2', 'Weighted_Mean2', 'p_value2'])
        for snp in beta_hash.keys():
            beta_list = beta_hash[snp]
            se_list = se_hash[snp]
            beta_list2 = beta_hash2[snp]
            se_list2 = se_hash2[snp]
            chrom_i = chrom_hash[snp][0]
            pos_i = pos_hash[snp][0]

            result = [snp, chrom_i, pos_i] + list(altmeta(beta_list, se_list)) + list(altmeta(beta_list2, se_list2))
            result_df.loc[len(result_df)] = result

        return result_df

    else:
        for _, row in file.iterrows():
            snp = row['SNP']
            beta, se = row[3], row[4]
            chrom, pos = row[1], row[2]
            beta_hash.setdefault(snp, []).append(beta)
            se_hash.setdefault(snp, []).append(se)
            chrom_hash.setdefault(snp, []).append(chrom)
            pos_hash.setdefault(snp, []).append(pos)

        result_df = pd.DataFrame(columns=['SNP', 'CHR', 'BP', 'P', 'I2', 'tau2_DL', 'p_Q', 'se', 'Z', 'Weighted_Mean', 'Q'])

        for snp in beta_hash.keys():
            beta_list = beta_hash[snp]
            se_list = se_hash[snp]
            chrom_i = chrom_hash[snp][0]
            pos_i = pos_hash[snp][0]

            result = [snp, chrom_i, pos_i] + list(altmeta(beta_list, se_list))
            result_df.loc[len(result_df)] = result

        return result_df
