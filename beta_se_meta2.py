import dask.dataframe as dd
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats.distributions import chi2


def altmeta(y1, s2):
    """Performs meta-analysis calculations."""
    s2 = np.array(s2, dtype=np.float64)
    variances = s2 ** 2
    weights = np.divide(1, s2, where=s2 > 0)  # Avoid division by zero
    n = len(y1)

    if np.isclose(weights.sum(), 0):
        return (None, 0, 0, None, None, None, None, 0)

    # Weighted mean and Q statistic
    mu_bar = np.dot(weights, y1) / weights.sum()
    Q = np.dot(weights, (y1 - mu_bar) ** 2)
    p_Q = chi2.sf(Q, n - 1)

    # Heterogeneity metrics
    I2 = max(0, (Q - (n - 1)) / Q) if Q > 0 else 0
    denom = weights.sum() - (weights ** 2).sum() / weights.sum()
    tau2_DL = max(0, (Q - n + 1) / denom) if denom > 0 else 0

    # Adjusted weights
    adjusted_variances = variances + tau2_DL
    adjusted_weights = np.divide(1, adjusted_variances, where=adjusted_variances > 0)

    if np.isclose(adjusted_weights.sum(), 0):
        return (None, I2, tau2_DL, p_Q, None, None, None, Q)

    # Final calculations
    mu_bar_adjusted = np.dot(adjusted_weights, y1) / adjusted_weights.sum()
    se = np.sqrt(1 / adjusted_weights.sum())
    z = mu_bar_adjusted / se if se > 0 else 0
    p = norm.sf(abs(z)) * 2

    return (p, I2, tau2_DL, p_Q, se, z, mu_bar_adjusted, Q)
  
  
# Define a helper function for processing each SNP group
def process_group(group):
    beta_list = group['BETA'].tolist()
    se_list = group['SE'].tolist()
    chrom = group['CHR'].iloc[0]
    pos = group['BP'].iloc[0]
    
    # Perform meta-analysis (placeholder for altmeta)
    altmeta_result = altmeta(beta_list, se_list)
    
    return [group.name, chrom, pos] + list(altmeta_result)

def beta_se(file):
    """Processes single beta and SE data."""
     # Group by SNP to process related data together
    grouped = file.groupby('SNP')

   

    # Apply the helper function to each group
    # Apply the helper function to each group
    results = grouped.apply(lambda g: pd.Series(process_group(g)))

    # Convert to DataFrame
    columns = ['SNP', 'CHR', 'BP', 'P', 'I2', 'tau2_DL', 'p_Q', 'SE', 'Z', 'Weighted_Mean', 'Q']
    result_df = pd.DataFrame(results, columns=columns)

    # Convert to Dask DataFrame with optimized partitioning
    return dd.from_pandas(result_df, npartitions=8)

def biv_beta_se(file):
    """Processes two betas and SEs data."""
    beta1_hash, se1_hash, beta2_hash, se2_hash, chrom_hash, pos_hash = {}, {}, {}, {}, {}, {}

    # Populate hash tables
    for row in file.itertuples(index=False):
        snp = row.SNP
        beta1, se1, beta2, se2 = row.BETA1, row.SE1, row.BETA2, row.SE2
        chrom, pos = row.CHR, row.BP

        beta1_hash.setdefault(snp, []).append(beta1)
        se1_hash.setdefault(snp, []).append(se1)
        beta2_hash.setdefault(snp, []).append(beta2)
        se2_hash.setdefault(snp, []).append(se2)
        chrom_hash[snp] = chrom
        pos_hash[snp] = pos

    # Process results
    results = []
    for snp, beta1_list in beta1_hash.items():
        se1_list = se1_hash[snp]
        beta2_list = beta2_hash[snp]
        se2_list = se2_hash[snp]
        chrom = chrom_hash[snp]
        pos = pos_hash[snp]

        # Perform meta-analysis for both sets of betas and SEs
        altmeta_result1 = altmeta(beta1_list, se1_list)
        altmeta_result2 = altmeta(beta2_list, se2_list)
        results.append((snp, chrom, pos) + altmeta_result1 + altmeta_result2)

    # Define columns
    columns = ['SNP', 'CHR', 'BP',
               'P1', 'I2_1', 'tau2_DL_1', 'p_Q1', 'SE1', 'Z1', 'Weighted_Mean1', 'Q1',
               'P2', 'I2_2', 'tau2_DL_2', 'p_Q2', 'SE2', 'Z2', 'Weighted_Mean2', 'Q2']

    # Create Dask DataFrame
    result_df = pd.DataFrame(results, columns=columns)
    return dd.from_pandas(result_df, npartitions=8)
