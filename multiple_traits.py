import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np
import math
import argparse
from pathlib import Path
from scipy.stats import t
from scipy.stats import chi2
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import pearsonr
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.special import chdtrc as chi2_cdf
from scipy.stats import moyal
from collections import defaultdict

 
def multiple_traits(data_dir,output_dir,method,keep_all ):
  def logit(p_values):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all probabilities are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
  
      k = len(p_values)
      C = np.sqrt(k * np.pi ** 2 * (5 * k + 2) / (3 * (5 * k + 4)))
      # Compute the logit transform
      t_value = -np.sum(np.log((p_values) / (1 - p_values))) / C
      df = 2 * k
      combined_p = 2 * (1 - t.cdf(np.abs(t_value), df))
  
      return combined_p
  
  
  def meanp(p_values):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all p-values are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-09, 1 - 1e-09)
      k = len(p_values)
      # Compute the mean of the p-values
      z_value = (0.5 - np.mean(p_values)) * math.sqrt(12 * k)
      combined_p = 1 - norm.cdf(z_value)  # right-side test
      return combined_p
  
  
  def fisher_method(p_values):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all p-values are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
      # Compute the combined test statistic using Fisher's method
      chi_squared = -2 * np.sum(np.log(p_values))
      # Calculate the degrees of freedom
      df = 2 * len(p_values)
  
      # Calculate the combined p-value using the chi-squared distribution
      fisher_p = 1 - chi2.cdf(chi_squared, df)
  
      return fisher_p
  
  
  def stouffer(p_values, one_tailed=True):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all p-values are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
      if one_tailed:
          z_scores = norm.ppf(1 - p_values)  # norm.ppf =inverse cumulative distribution function of normal distribution
      else:
          z_scores = norm.ppf(1 - p_values / 2)
  
      # z_scores = norm.ppf(1 - p_values) # norm.ppf =inverse cumulative distribution function of normal distribution
      combined_z = np.sum(z_scores / math.sqrt(len(p_values)))
      combined_p = 1 - norm.cdf(combined_z)  # norm.cdf = cumulative distribution function of normal distribution
  
      return combined_p
  
  
  def inverse_chi2(p_values):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all p-values are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
  
      df = len(p_values)
  
      chi_squared = np.sum(chi2.ppf((1 - p_values), 1))  # ppf = inverse of cdf
      inv_chi2_combined_p = 1 - chi2.cdf(chi_squared, df)  # chi2.cdf = cumulative distribution function of chi2
  
      return inv_chi2_combined_p
  
  
  def binomial_test(p_values):
      p_values = np.array(p_values)
  
      k = len(p_values)
      alpha = 0.05
      # Count the number of significant p-values
      r = sum((p < alpha) for p in p_values)
      # Calculate the binomial probability of observing at most num_successes successes
      combined_p = 0
      for x in range(r, k + 1):
          combined_p += math.factorial(k) / (math.factorial(x) * math.factorial(k - x)) * (alpha ** x) * (
                      (1 - alpha) ** (k - x))
  
      return combined_p
  
  
  def cauchy_cdf(x, x0=0, gamma=1):
      """Cumulative distribution function (CDF) of the Cauchy distribution."""
  
      return 0.5 + (np.arctan((x - x0) / gamma) / np.pi)
  
  
  def cauchy_method(p_values):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all p-values are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-15, 1 - 1e-15)
  
      k = len(p_values)
  
      T = np.tan((0.5 - p_values) * np.pi)
  
      t = np.sum(T) / k
      # Calculate the combined p-value using the Cauchy distribution
      combined_p = 1 - cauchy_cdf(t)
  
      return combined_p
  
  
  def minP(p_values):
      # Convert input to numpy array for numerical operations
      p_values = np.array(p_values)
  
      # Ensure that all p-values are within the valid range (0, 1)
      p_values = np.clip(p_values, 1e-15, 1 - 1e-15)
  
      k = len(p_values)
  
      # Calculate the combined p-value as the minimum p-value
      combined_p = 1 - (1 - np.min(p_values)) ** k
  
      return combined_p
  
  
  def CMC(p_values):
      p_value_cauchy = cauchy_method(p_values)
      p_value_minp = minP(p_values)
      combined_p = cauchy_method((p_value_cauchy, p_value_minp))  # pCMC = CCT{pCCT , pMinP}
      return combined_p
  
  
  def MCM(p_values):
      p_value_cauchy = cauchy_method(p_values)
      p_value_minp = minP(p_values)
      combined_p = 2 * min(p_value_cauchy, p_value_minp, 0.5)  # pMCM = 2 min{pCCT , pMinP, 0.5}
  
      return combined_p
  
  
  def HMP(p_values):
      p_values = np.array(p_values)
      L = len(p_values)
      harmonic_mean = L / np.sum(1 / p_values)
      combined_p = moyal.cdf(-2 * np.log(harmonic_mean), 1, L)
      return combined_p
  
  
  # Function to combine p-values based on specified method
  def combine_p_values(p_values, method):
      if method == 'logit':
          return logit(p_values)
      elif method == 'meanp':
          return meanp(p_values)
      elif method == 'fisher':
          return fisher_method(p_values)
      elif method == 'stouffer':
          return stouffer(p_values)
      elif method == 'inverse_chi2':
          return inverse_chi2(p_values)
      elif method == 'binomial_test':
          return binomial_test(p_values)
      elif method == 'cauchy':
          return cauchy_method(p_values)
      elif method == 'minP':
          return minP(p_values)
      elif method == 'CMC':
          return CMC(p_values)
      elif method == 'MCM':
          return MCM(p_values)
      elif method == 'HMP':
          return HMP(p_values)
      else:
          raise ValueError(f"Unknown method: {method}")
  
  significance_threshold = 1e-7
  # Initialize list for merging p-values
  all_study_data = []
  
  # Loop through each study file
  for file_path in data_dir.glob("*.txt"):
      # Read each study file
      study_df = pd.read_csv(file_path, sep='\t')
      
      ############ Imputation part here ################
      
      
      
      ############ Imputation part here ################
      # Extract filename for naming
      filename = file_path.stem
      
      # Precompute -log10(P) values with minimum cap to avoid -inf
      study_df['-log10(P)'] = -np.log10(study_df['P'].clip(lower=1e-300))
      
      # Plot Manhattan plot
      plt.figure(figsize=(15, 6))
      
      # Define colors for chromosomes
      unique_chromosomes = study_df['CHR'].unique()
      colors = ['black', 'grey'] * (len(unique_chromosomes) // 2 + 1)
      
      # Calculate positions and xticks for each chromosome
      chrom_positions = study_df.groupby('CHR')['BP'].max().cumsum() + np.arange(len(unique_chromosomes)) * 1e6
      chrom_start_positions = dict(zip(unique_chromosomes, chrom_positions.shift(fill_value=0)))
      
      # Scatter plot chromosomes with precomputed chromosomal positions
      for chrom, color in zip(unique_chromosomes, colors):
          subset = study_df[study_df['CHR'] == chrom]
          plt.scatter(subset['BP'] + chrom_start_positions[chrom], subset['-log10(P)'], color=color, s=1)
      
      # Threshold line and plot styling
      plt.axhline(y=-np.log10(significance_threshold), color='red', linestyle='--')
      plt.xlabel('Chromosome')
      plt.ylabel('-log10(P)')
      plt.title(f'Manhattan Plot for {filename}')
      plt.xticks(ticks=chrom_positions, labels=unique_chromosomes)
      plt.savefig(f"plots/{filename}_manhattan.png", dpi=600)
      plt.close()
      
      # Append study data for final merging
      all_study_data.append(study_df[['SNP', 'CHR', 'BP', 'P']].rename(columns={'P': f'P_val_from_{filename}'}))
  
  # Concatenate all studies' data and merge
  # Concatenate all studies' data and merge based on SNP
  final_combined_df = all_study_data[0]
  for additional_df in all_study_data[1:]:
    final_combined_df = final_combined_df.merge(additional_df, on=['SNP', 'CHR', 'BP'], how='outer')
  
  if keep_all == 'NO':
    final_combined_df   = final_combined_df.dropna()
 
      
  print( final_combined_df.head())
  
  # Extract the necessary columns starting from the third column onward (p-values)
  values = np.array(final_combined_df.iloc[:, 3:].values)
  
  # Create a copy of `final_combined_df` without the individual p-value columns
  final_combined_df_no_pvalues = final_combined_df.iloc[:, :3].copy()  # Retain only the first three columns
  
  # Apply the selected combination method to each row of p-values and add it as a new 'P' column
  final_combined_df_no_pvalues['P'] = np.apply_along_axis(combine_p_values, 1, values, method=method)
  # Save the final DataFrame without individual p-value columns to CSV
  final_combined_df_no_pvalues.to_csv(output_dir / f"multiple_traits_results.txt", index=False, sep='\t')
  
  #################### Manhattan plot here for the combined p #########################
  # Precompute -log10(P) values with minimum cap to avoid -inf
  final_combined_df_no_pvalues['-log10(P)'] = -np.log10(final_combined_df_no_pvalues['P'].clip(lower=1e-300))
  
  # Plot Manhattan plot
  plt.figure(figsize=(15, 6))
  
  # Define colors for chromosomes
  unique_chromosomes = final_combined_df_no_pvalues['CHR'].unique()
  colors = ['black', 'grey'] * (len(unique_chromosomes) // 2 + 1)
  
  # Calculate positions and xticks for each chromosome
  chrom_positions = final_combined_df_no_pvalues.groupby('CHR')['BP'].max().cumsum() + np.arange(len(unique_chromosomes)) * 1e6
  chrom_start_positions = dict(zip(unique_chromosomes, chrom_positions.shift(fill_value=0)))
  
  # Scatter plot chromosomes with precomputed chromosomal positions
  for chrom, color in zip(unique_chromosomes, colors):
      subset = final_combined_df_no_pvalues[final_combined_df_no_pvalues['CHR'] == chrom]
      plt.scatter(subset['BP'] + chrom_start_positions[chrom], subset['-log10(P)'], color=color, s=1)
  
  # Threshold line and plot styling
  plt.axhline(y=-np.log10(significance_threshold), color='red', linestyle='--')
  plt.xlabel('Chromosome')
  plt.ylabel('-log10(P)')
  plt.title(f'Manhattan Plot - Combined p-value {method}')
  plt.xticks(ticks=chrom_positions, labels=unique_chromosomes)
 
  
  # Save the plot with 600 DPI
  plt.savefig(f"plots/manhattan_plot_combined.png", dpi=600)
  plt.close()
  #################### Manhattan plot here for the combined p #########################

  
  


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine p-values and generate Manhattan plots for multiple traits.")
    parser.add_argument("--data_dir", type=str, default="studies", help="Directory containing study files.")
    parser.add_argument("--output_dir", type=str, default="plots", help="Directory to save output plots and combined p-values.")
    parser.add_argument("--method", type=str,default="fisher",help=(
                              "Method for combining p-values. Available options are:\n"
                              "  - fisher: Fisher's method\n"
                              "  - logit: Logit method\n"
                              "  - meanp: Mean of p-values method\n"
                              "  - stouffer: Stouffer's method (one-tailed)\n"
                              "  - inverse_chi2: Inverse chi-squared method\n"
                              "  - binomial_test: Binomial test method\n"
                              "  - cauchy: Cauchy combination method\n"
                              "  - minP: Minimum p-value method\n"
                              "  - CMC: Combined Cauchy and minP (pCMC) method\n"
                              "  - MCM: Minimum Cauchy and minP method\n"
                              "  - HMP: Harmonic mean p-value method"
                          )
                      )
                      
    parser.add_argument("--keep_all", type=str, default="NO",   help="Specify whether to retain all variants, even those missing from the study. Set to 'YES' to keep all variants, or 'NO' (default) to filter out missing ones.")

    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    multiple_traits(data_dir, output_dir, args.method,args.keep_all)

 
