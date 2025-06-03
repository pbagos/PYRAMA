import os
import subprocess
import pandas as pd
import shutil
import sys
import dask.dataframe as dd 
def imputation(data_file, output_folder, r2threshold, population, maf, ref, imp_list=None):
    print("Running Imputation:")
    os.makedirs(output_folder, exist_ok=True)

    command = [
        "python3", "pred_ld.py",
        "--file-path", data_file,
        "--r2threshold", str(r2threshold),
        "--pop", str(population),
        "--maf", str(maf),
        "--ref", str(ref)
    ]
    if imp_list:
        command += ["--imp_list", imp_list]

    try:
        subprocess.run(command, check=True)
        print("Imputation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during imputation: {e}")
        return None

    # Move files to output folder
    for filename in os.listdir():
        if filename.startswith("LD_info") or filename.startswith("imputation_results"):
            shutil.move(filename, os.path.join(output_folder, filename))

    # Load results into DataFrame
    imputation_results_files = [f for f in os.listdir(output_folder) if f.startswith("imputation_results")]
    if imputation_results_files:
        imputation_results_path = os.path.join(output_folder, imputation_results_files[0])
        data = pd.read_csv(imputation_results_path, sep='\t')
        
        output_folder = "results"
        # Step 1: Identify SNPs that exist both as imputed (1) and non-imputed (0)
        duplicate_snps = data[data['SNP'].duplicated(keep=False)]['SNP'].unique()
        
        # Step 2: Filter the data
        # Keep non-imputed rows for duplicate SNPs and keep all other rows as they are
        filtered_data = data[~data['SNP'].isin(duplicate_snps) | (data['imputed'] == 0)]
        
        # Display or save the filtered data
        data = filtered_data 
        
        data = data[['SNP','CHR','BP','BETA','SE']]
        
        # Ensuring correct data types
        data['SNP'] = data['SNP'].astype(str)          # SNP as string
        data['CHR'] = data['CHR'].astype(int)          # CHR as integer
        #data['BP'] = data['BP'].astype(int)            # BP as integer
        data['BETA'] = data['BETA'].astype(float)      # BETA as float
        data['SE'] = data['SE'].astype(float)          # SE as float
        data = data.dropna() 
        data.to_csv("results/imputation_results.txt", sep='\t', index=False)
        print(f"Imputation results saved to results/imputation_results.txt")
        return data
    else:
        print("No imputation results found.")
        return None



if __name__ == "__main__":
    imputation(*sys.argv[1:])
