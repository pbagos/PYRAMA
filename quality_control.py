# Imports
import pandas as pd
import argparse
import os
import warnings
from itertools import combinations

# Settings the warnings to be ignored
warnings.filterwarnings('ignore')


def read_gwas_file(filepath):
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

    expected_cols = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE']
    mandatory_cols = ['SNP', 'BETA', 'SE']
    
    missing_cols = [col for col in expected_cols if col not in df.columns]
    missing_cols_mandatory = [col for col in mandatory_cols if col not in df.columns]

    if missing_cols:
        print(f"File {filepath} is missing columns: {missing_cols}")
    if missing_cols_mandatory:
        print(f"File {filepath} is missing mandatory columns: {missing_cols_mandatory}")
        return None

    present_expected_cols = [col for col in expected_cols if col in df.columns]
    df = df[present_expected_cols]

    return df


def quality_control(df, skip_harm=False):
    problematic_lines = []
    initial_count = len(df)

    missing_mask = df[['BETA', 'SE']].isnull().any(axis=1)
    missing_rows = df[missing_mask]
    for idx, row in missing_rows.iterrows():
        problematic_lines.append({
            "study": row["study"],
            "SNP": row["SNP"],
            "reason": "Missing BETA or SE"
        })
    df = df.dropna(subset=['BETA', 'SE'])

    df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')
    df['SE'] = pd.to_numeric(df['SE'], errors='coerce')

    conversion_mask = df[['BETA', 'SE']].isnull().any(axis=1)
    conversion_rows = df[conversion_mask]
    for idx, row in conversion_rows.iterrows():
        problematic_lines.append({
            "study": row["study"],
            "SNP": row["SNP"],
            "reason": "Invalid numeric value in BETA or SE"
        })
    df = df.dropna(subset=['BETA', 'SE'])

    negative_se_mask = df['SE'] <= 0
    negative_se_rows = df[negative_se_mask]
    for idx, row in negative_se_rows.iterrows():
        problematic_lines.append({
            "study": row["study"],
            "SNP": row["SNP"],
            "reason": "Wrong SE value"
        })
    df = df[~negative_se_mask]

    if not skip_harm and 'A1' in df.columns and 'A2' in df.columns:
        valid_alleles = {"A", "G", "C", "T"}
        allele_mask = df['A1'].isin(valid_alleles) & df['A2'].isin(valid_alleles)
        allele_problematic_rows = df[~allele_mask]
        for idx, row in allele_problematic_rows.iterrows():
            problematic_lines.append({
                "study": row["study"],
                "SNP": row["SNP"],
                "reason": "wrong allele symbol"
            })
        df = df[allele_mask]

    final_count = len(df)
    print(f"Quality Control: Dropped {initial_count - final_count} rows due to missing/invalid BETA, SE values or wrong allele symbol")

    prob_df = pd.DataFrame(problematic_lines)
    if not prob_df.empty:
        prob_df = prob_df[['study', 'SNP', 'reason']]
    else:
        prob_df = pd.DataFrame(columns=['study', 'SNP', 'reason'])

    return df, initial_count, final_count, prob_df


def harmonize_variants(df):
    harmonization_count = 0
    grouped = df.groupby('SNP')
    for snp, group in grouped:
        ref_allele = group.iloc[0]['A1']
        for idx in group.index:
            current_a1 = df.at[idx, 'A1']
            if current_a1 != ref_allele:
                df.at[idx, 'BETA'] = -df.at[idx, 'BETA']
                original_a1 = df.at[idx, 'A1']
                df.at[idx, 'A1'] = df.at[idx, 'A2']
                df.at[idx, 'A2'] = original_a1
                harmonization_count += 1
    return df, harmonization_count


def compute_common_snps(df):
    unique_studies = df["study"].unique()
    for study in unique_studies:
        print(f"SNPs that exist in {study} and not in the other studies are written to {study}_missing_snps_from_all_others.txt")

    study_combinations = []
    for r in range(1, len(unique_studies) + 1):
        study_combinations.extend(combinations(unique_studies, r))

    results = []
    for combo in study_combinations:
        row = {study: "X" if study in combo else "O" for study in unique_studies}
        common_snps = set(df[df["study"] == combo[0]]["SNP"])
        for study in combo[1:]:
            common_snps.intersection_update(set(df[df["study"] == study]["SNP"]))
        row["number_of_SNPs"] = len(common_snps)
        results.append(row)

        for study in unique_studies:
            snps_study = set(df[df["study"] == study]["SNP"])
            other_studies = [s for s in unique_studies if s != study]
            snps_all_others = set().union(*(set(df[df["study"] == s]["SNP"]) for s in other_studies))
            missing_snps_combined = snps_all_others - snps_study

            filename_combined = f"{study}_missing_snps_from_all_others.txt"
            with open(filename_combined, "w") as file:
                file.write("\n".join(missing_snps_combined))

    return pd.DataFrame(results)


def filter_consistent_snps(df):
    df['allele_pair'] = df.apply(lambda row: tuple(sorted([row['A1'], row['A2']])), axis=1)
    def is_consistent(group): return group['allele_pair'].nunique() == 1
    consistent_df = df.groupby('SNP').filter(is_consistent)
    inconsistent_df = df.groupby('SNP').filter(lambda group: not is_consistent(group))
    consistent_df = consistent_df.drop('allele_pair', axis=1)
    problematic_df = inconsistent_df[['study', 'SNP']].copy()
    problematic_df['reason'] = 'inconsistent alleles'
    return consistent_df, problematic_df


def main():
    print("-----------------------------------------------------")
    print("Quality Control")
    print("-----------------------------------------------------")

    parser = argparse.ArgumentParser(
        description="Combine multiple GWAS files, perform quality control, and allele harmonization.")
    parser.add_argument('input_files', nargs='+',
                        help='Input GWAS files (tab-delimited) with columns: SNP, CHR, BP, A1, A2, BETA, SE.')
    parser.add_argument('--output', default='merged_study.txt',
                        help='Output file for merged and harmonized GWAS data.')
    parser.add_argument('--skip_harm', action='store_true',
                        help='Skip allele harmonization and ignore A1/A2 filtering.')
    args = parser.parse_args()

    snp_counts = {}
    dataframes = []
    problematic_dfs = []

    for filepath in args.input_files:
        print(f"Processing file: {filepath}")
        df = read_gwas_file(filepath)
        if df is None:
            print(f"Skipping file: {filepath} due to errors.")
            continue
        study_name = os.path.basename(filepath)
        df['study'] = study_name
        unique_snps = df['SNP'].nunique()
        snp_counts[study_name] = unique_snps
        dataframes.append(df)

    if not dataframes:
        print("No valid GWAS files provided. Exiting.")
        return

    merged_df = pd.concat(dataframes, ignore_index=True)
    print(f"Total merged records before quality control: {len(merged_df)}")

    merged_df, initial_count, final_count, prob_df = quality_control(merged_df, skip_harm=args.skip_harm)
    merged_df = merged_df.sort_values(by=['SNP', 'CHR'])

    harmonization_count = 0
    problematic_df_allele = pd.DataFrame(columns=['study', 'SNP', 'reason'])

    if not args.skip_harm and "A1" in merged_df.columns and "A2" in merged_df.columns:
        merged_df, problematic_df_allele = filter_consistent_snps(merged_df)
        merged_df, harmonization_count = harmonize_variants(merged_df)
        print(f"Total allele harmonizations performed: {harmonization_count}")
    elif args.skip_harm:
        print("Harmonization skipped as requested.")
    else:
        print("Allele columns not found. Harmonization not performed.")

    print(f"Total merged records after quality control: {len(merged_df)}")

    for study in merged_df['study'].unique():
        df_study = merged_df[merged_df['study'] == study]
        output_file = study.split(".")[0] + "_harmonized.txt"
        df_study = df_study.drop('study', axis=1)
        df_study.to_csv(output_file, sep='\t', index=False)
        print(f"Written {len(df_study)} rows to {output_file}")

    merged_df.to_csv(args.output, sep='\t', index=False)
    print(f"Merged and harmonized GWAS data written to {args.output}")

    prob_df = pd.concat([prob_df, problematic_df_allele])
    if not prob_df.empty:
        prob_df.to_csv("prob_lines.txt", sep='\t', index=False)
        print("Problematic lines written to prob_lines.txt")

    common_array = compute_common_snps(merged_df)
    common_array.to_csv("common_SNPs_combinations.txt", sep='\t', index=False)
    print("Common SNPs across all study combinations written to common_SNPs_combinations.txt")

    report_lines = []
    report_lines.append("Quality Control Report")
    report_lines.append("-----------------------\n")
    report_lines.append("Unique SNP counts per study:")
    for filename, count in snp_counts.items():
        report_lines.append(f"{filename}: {count}")
    report_lines.append(f"\nTotal allele harmonizations performed: {harmonization_count}")
    if not prob_df.empty:
        report_lines.append(f"\nDropped Lines: {initial_count - final_count}. See prob_lines.txt for further details\n")

    with open("quality_control_report.txt", "w") as report_file:
        report_file.write("\n".join(report_lines))
    print("Quality control report written to quality_control_report.txt")
    print("-----------------------------------------------------")


if __name__ == '__main__':
    main()
