import os
import pandas as pd
import numpy as np

import cont_model
import argparse
import polars as pl
from scipy.stats import norm, cauchy, pearsonr
from decimal import Decimal, getcontext

# Set the precision for Decimal calculations
getcontext().prec = 1000


def cauchy_cdf(x, x0=0, gamma=1):
    """Cumulative distribution function (CDF) of the Cauchy distribution."""
    return 0.5 + (np.arctan((x - x0) / gamma) / np.pi)


def random_effects_meta_analysis(effects, variances, alpha=0.05, het_est="DL"):
    effects = np.array(effects)
    variances = np.array(variances)
    heterogeneity = het_est

    weights_fixed = 1.0 / variances
    fixed_effect = np.sum(weights_fixed * effects) / np.sum(weights_fixed)
    Q = np.sum(weights_fixed * (effects - fixed_effect) ** 2)
    k = len(effects)
    df = k - 1
    denom = np.sum(weights_fixed) - np.sum(weights_fixed ** 2) / np.sum(weights_fixed)

    # Heterogeneity estimator selection
    if heterogeneity.upper() == "DL":
        tau2 = max(0, (Q - df) / denom) if denom > 0 else 0.0

    elif heterogeneity.upper() == "ANOVA":
        numerator = np.sum((effects - np.mean(effects)) ** 2) / (k - 1)
        within_variance_avg = np.sum(variances) / k
        tau2 = max(0, numerator - within_variance_avg)

    elif heterogeneity.upper() == "SJ":

        yi = effects
        vi = variances

        # Number of studies
        k = yi.size

        X = np.ones((k, 1))
        p = X.shape[1]

        Y_bar = yi.mean()  # sample mean of betas
        ymci = yi - Y_bar  # mean centered effects

        tau2_0 = np.var(ymci, ddof=0)

        wi = 1.0 / (vi + tau2_0)
        W = np.diag(wi)  # (k × k)

        XtWX = X.T @ W @ X
        inv_XtWX = np.linalg.inv(XtWX)

        P = W - W @ X @ inv_XtWX @ X.T @ W

        beta_hat = inv_XtWX @ (X.T @ W @ yi)  # shape (1,)

        fitted = X @ beta_hat  # shape (k,)
        Ymc = yi - fitted  # residuals

        RSS = float(Ymc.T @ P @ Ymc)

        tau2 = tau2_0 * RSS / (k - p)

    elif heterogeneity.upper() == "FE":
        tau2 = 0.0

    else:
        raise ValueError(f"Unsupported heterogeneity estimator: {heterogeneity}. Use 'DL', 'ANOVA', 'SJ', or 'FE'.")

    weights_random = 1.0 / (variances + tau2)
    overall_effect = np.sum(weights_random * effects) / np.sum(weights_random)
    overall_se = np.sqrt(1.0 / np.sum(weights_random))
    return overall_effect / overall_se


def correlated_Stouffer(values, cor_sum, k):
    values = np.array(values)
    z_scores = values

    total_variance = k + 2 * cor_sum
    # Compute combined test statistic
    combined_z = np.sum(z_scores) / np.sqrt(total_variance)

    # Compute combined p-value
    combined_p_value = 1 - norm.cdf(combined_z)

    return combined_p_value


def fast_robust_analysis(data,   het_est):
    # Drop rows with any missing values
    data = data.dropna()

    # Extract CHR and BP mapping for each SNP (takes the first occurrence per SNP)
    chr_bp = (
        data[['SNP', 'CHR', 'BP']]
        .drop_duplicates(subset='SNP')
        .set_index('SNP')
    )

    # Prepare data for meta-analysis (exclude SNP column)
    file = data.copy()
    file.index = file['SNP']
    file = file.drop(['SNP', 'CHR', 'BP'], axis=1)

    snp_keys = file.index.unique()
    snp_hash = {snp: file.loc[snp].values.tolist() for snp in snp_keys}

    # Containers for results
    snps_ = []
    z_dom_list, z_add_list, z_rec_list = [], [], []
    p_dom, p_add, p_rec = [], [], []
    p_value_min_p, p_value_cauchy = [], []

    # Loop through each SNP
    for snp_name in snp_keys:
        dom_meta_es, add_meta_es, rec_meta_es = [], [], []
        dom_meta_var, add_meta_var, rec_meta_var = [], [], []

        # Build a small DataFrame per SNP
        snp_data = pl.DataFrame(snp_hash[snp_name]).transpose()

        # Iterate over rows (each row represents one study)
        for row in snp_data.iter_rows():


            XAA,SDAA,NAA, XAB,SDAB,NAB, XBB,SDBB,NBB = [np.array(row[col]) for col in range(0, 9)]



            row_list = [ XAA,SDAA,NAA, XAB,SDAB,NAB, XBB,SDBB,NBB]

            # Compute effect sizes and variances
            effect_dom, var_dom = cont_model.cont_model_dominant(row_list)
            effect_add, var_add = cont_model.cont_model_additive(row_list)
            effect_rec, var_rec = cont_model.cont_model_recessive(row_list )

            dom_meta_es.append(effect_dom)
            dom_meta_var.append(var_dom)
            add_meta_es.append(effect_add)
            add_meta_var.append(var_add)
            rec_meta_es.append(effect_rec)
            rec_meta_var.append(var_rec)

        # Random-effects meta-analysis for each genetic model
        z_dom = random_effects_meta_analysis(dom_meta_es, dom_meta_var, alpha=0.05, het_est=het_est)
        z_add = random_effects_meta_analysis(add_meta_es, add_meta_var, alpha=0.05, het_est=het_est)
        z_rec = random_effects_meta_analysis(rec_meta_es, rec_meta_var, alpha=0.05, het_est=het_est)

        z_dom_list.append(z_dom)
        z_add_list.append(z_add)
        z_rec_list.append(z_rec)

        # Individual p-values
        p_vals = np.array([norm.sf(abs(z_dom)) * 2,
                           norm.sf(abs(z_add)) * 2,
                           norm.sf(abs(z_rec)) * 2])

        # MinP combination
        p_min = p_vals.min()
        combined_p = 1 - (1 - Decimal(p_min)) ** 3
        p_value_min_p.append(np.float64(combined_p))

        # Cauchy combination
        T = np.tan((0.5 - p_vals) * np.pi)
        t = T.mean()
        p_cauchy = Decimal(cauchy.sf(t))
        p_value_cauchy.append(np.float64(p_cauchy))

        # Store per-model p-values
        p_dom.append(p_vals[0])
        p_add.append(p_vals[1])
        p_rec.append(p_vals[2])
        snps_.append(snp_name)

    # Convert combined p-lists to arrays
    p_value_min_p = np.array(p_value_min_p)
    p_value_cauchy = np.array(p_value_cauchy)

    # CMC (Combined Cauchy-MinP)
    tan_terms = np.tan((0.5 - p_value_cauchy) * np.pi) + np.tan((0.5 - p_value_min_p) * np.pi)
    p_cmc = 0.5 - np.arctan(tan_terms / 2) / np.pi

    # MCM (Minimum Cauchy-MinP)
    min_vals = np.minimum(p_value_cauchy, p_value_min_p)
    p_mcm = np.minimum(1, 2 * min_vals)

    # Correlated Stouffer
    stouffer_corr_p = []
    data_matrix = np.column_stack((z_dom_list, z_add_list, z_rec_list))

    k = data_matrix.shape[1]  # Number of tests (columns)
    # Number of SNPs (observations)
    n_snps = len(z_dom_list)
    # Prepare storage for the Stouffer p-values
    if n_snps <= 1:
        # If only one SNP, we cannot estimate correlations – fill with NaN
        stouffer_corr_p = [np.nan] * n_snps
    else:
        # Stack your Z‐lists into an (n_snps × 3) matrix
        data_matrix = np.column_stack((z_dom_list, z_add_list, z_rec_list))
        k = data_matrix.shape[1]  # should be 3 tests: Dom, Add, Rec

        # Compute sum of pairwise correlations across tests
        cor_sum = 0.0
        for i in range(k):
            for j in range(i + 1, k):
                cor, _ = pearsonr(data_matrix[:, i], data_matrix[:, j])
                cor_sum += cor

        # Apply correlated Stouffer to each row
        stouffer_corr_p = []
        for row in data_matrix:
            p_val = correlated_Stouffer(values=row, cor_sum=cor_sum, k=k)
            stouffer_corr_p.append(p_val)

    # Build result DataFrame
    result = pd.DataFrame({
        'SNP':       snps_,
        'Z_Dom':     z_dom_list,
        'Z_Add':     z_add_list,
        'Z_Rec':     z_rec_list,
        'P_Dom':     np.array(p_dom),
        'P_Add':     np.array(p_add),
        'P_Rec':     np.array(p_rec),
        'P_MinP':    p_value_min_p,
        'P_CCT':     p_value_cauchy,
        'P_CMC':     p_cmc,
        'P_MCM':     p_mcm,
        'P_Stouffer': stouffer_corr_p
    })

    # Merge in CHR/BP and reorder
    result = result.merge(chr_bp, left_on='SNP', right_index=True)
    cols = ['SNP', 'CHR', 'BP'] + [c for c in result.columns if c not in ['SNP', 'CHR', 'BP']]
    result = result[cols]

    return result



if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Fast robust analysis for genetic data.")

    parser.add_argument("--path", type=str, required=True, help="Path to the directory containing input files.")
    parser.add_argument("--file_list", type=str, required=True, nargs='+',
                        help="List of input data files (space-separated).")
    parser.add_argument("--output", type=str, default="results.csv", help="Output file name for saving the results.")

    parser.add_argument('--het_est', required=False, default="DL", choices=["DL", "ANOVA", "SJ", "FE"],
                        help="Heterogeneity estimator: 'DL' (DerSimonian-Laird), 'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman), or 'FE' (Fixed Effects Only)")
    # Parse arguments
    args = parser.parse_args()

    # Perform the analysis
    final_results = fast_robust_analysis(path=args.path, file_list=args.file_list,
                                         effect_size_type=args.effect_size_type, het_est=args.het_est)

    # Save results to the output file
    final_results.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to {args.output}")

