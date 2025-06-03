# merge_studies.py

import os
import sys
import pandas as pd

def merge_studies(folder_path, output_combined, output_variants):
    print(f"Merging studies of the folder_path : {folder_path}")
    # List all study files in the specified folder path
    study_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".txt")]
    dfs = [pd.read_csv(f, delim_whitespace=True) for f in study_files]

    # Merge files on the SNP column and drop duplicates
    combined_df = pd.concat(dfs,ignore_index=True) 
    combined_df = combined_df.sort_values("SNP")
    print("Studies merged")
    print("Saving the merged_file")

    # Save combined data
    combined_df.to_csv(output_combined, sep="\t", index=False)

    # Find unique SNPs not common to all studies
    all_snps = set(combined_df["SNP"])
    non_common_snps = []
    for df in dfs:
        snps_in_study = set(df["SNP"])
        non_common_snps.extend(list(all_snps - snps_in_study))

    # Save non-common SNPs to the output file
    with open(output_variants, "w") as f:
        f.write("\n".join(non_common_snps))

if __name__ == "__main__":
    # Read arguments from command line
    folder_path = sys.argv[1]
    output_combined = sys.argv[2]
    output_variants = sys.argv[3]

    # Run the merge studies function
    merge_studies(folder_path, output_combined, output_variants)
