import argparse
import os
import warnings
import polars as pl
import pandas as pd 
from itertools import combinations
# Ignore warnings
warnings.filterwarnings('ignore')


def read_gwas_file(filepath: str) -> pl.DataFrame | None:
    """
    Read a GWAS file into a Polars DataFrame, checking for required columns.
    """
    try:
        df = pl.read_csv(filepath, separator="\t")
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

    expected_cols = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE']
    mandatory_cols = ['SNP', 'BETA', 'SE']

    missing_cols = [c for c in expected_cols if c not in df.columns]
    missing_mand = [c for c in mandatory_cols if c not in df.columns]

    if missing_cols:
        print(f"File {filepath} is missing columns: {missing_cols}")
    if missing_mand:
        print(f"File {filepath} is missing mandatory columns: {missing_mand}")
        return None

    # Select only expected columns present
    df = df.select([c for c in expected_cols if c in df.columns])
    return df


def quality_control(
    df: pl.DataFrame,
    skip_harm: bool = False
) -> tuple[pl.DataFrame, int, int, pl.DataFrame]:
    """
    Perform QC: remove missing/invalid BETA/SE, non-positive SE, and invalid alleles.
    Returns cleaned df, initial count, final count, and DataFrame of problematic lines.
    """
    problematic = []
    initial_count = df.height

    # 1) Missing BETA or SE
    missing = df.filter(pl.col('BETA').is_null() | pl.col('SE').is_null())
    for row in missing.select(['study', 'SNP']).to_dicts():
        problematic.append({**row, 'reason': 'Missing BETA or SE'})
    df = df.drop_nulls(['BETA', 'SE'])

    # 2) Ensure numeric
    df = df.with_columns([
        pl.col('BETA').cast(pl.Float64),
        pl.col('SE').cast(pl.Float64)
    ])
    invalid_numeric = df.filter(pl.col('BETA').is_null() | pl.col('SE').is_null())
    for row in invalid_numeric.select(['study', 'SNP']).to_dicts():
        problematic.append({**row, 'reason': 'Invalid numeric value in BETA or SE'})
    df = df.drop_nulls(['BETA', 'SE'])

    # 3) Non-positive SE
    nonpos_se = df.filter(pl.col('SE') <= 0)
    for row in nonpos_se.select(['study', 'SNP']).to_dicts():
        problematic.append({**row, 'reason': 'Wrong SE value'})
    df = df.filter(pl.col('SE') > 0)

    # 4) Allele filtering
    if not skip_harm and {'A1', 'A2'}.issubset(df.columns):
        valid = ['A', 'C', 'G', 'T']
        bad_allele = df.filter(~(pl.col('A1').is_in(valid) & pl.col('A2').is_in(valid)))
        for row in bad_allele.select(['study', 'SNP']).to_dicts():
            problematic.append({**row, 'reason': 'wrong allele symbol'})
        df = df.filter(pl.col('A1').is_in(valid) & pl.col('A2').is_in(valid))

    final_count = df.height
    print(
        f"Quality Control: Dropped {initial_count - final_count} rows "
        "due to missing/invalid BETA, SE values or wrong allele symbol"
    )

    prob_df = pl.DataFrame(problematic) if problematic else pl.DataFrame(
        {'study': [], 'SNP': [], 'reason': []}
    )
    return df, initial_count, final_count, prob_df


def harmonize_variants(
    df: pl.DataFrame
) -> tuple[pl.DataFrame, int]:
    """
    Harmonize allele signs so that A1 matches the reference from the first occurrence per SNP.
    """
    records = df.to_dicts()
    ref = {}
    harmon_count = 0
    # Determine reference A1 per SNP
    for rec in records:
        snp = rec['SNP']
        if snp not in ref:
            ref[snp] = rec['A1']
    # Flip BETA and swap alleles where needed
    for rec in records:
        if rec['A1'] != ref[rec['SNP']]:
            rec['BETA'] = -rec['BETA']
            rec['A1'], rec['A2'] = rec['A2'], rec['A1']
            harmon_count += 1
    return pl.DataFrame(records), harmon_count


def filter_consistent_snps(
    df: pl.DataFrame
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Keep only SNPs where allele pairs are consistent across studies.
    Return cleaned df and DataFrame of problematic SNPs.
    """
    # Create a canonical allele_pair string (sorted)
    df = df.with_columns(
        pl.when(pl.col('A1') < pl.col('A2'))
        .then(pl.concat_str([pl.col('A1'), pl.col('A2')], separator='_'))
        .otherwise(pl.concat_str([pl.col('A2'), pl.col('A1')], separator='_'))
        .alias('allele_pair')
    )
    # Determine SNPs with a single unique allele_pair
    counts = df.group_by('SNP').agg(
        pl.col('allele_pair').n_unique().alias('n_pairs')
    )
    good_snps = counts.filter(pl.col('n_pairs') == 1)['SNP'].to_list()

    consistent = df.filter(pl.col('SNP').is_in(good_snps)).drop('allele_pair')
    inconsistent = df.filter(~pl.col('SNP').is_in(good_snps))
    problem = inconsistent.select(['study', 'SNP']).with_columns(
        pl.lit('inconsistent alleles').alias('reason')
    )
    return consistent, problem


def compute_common_snps(
    df: pl.DataFrame
) -> pl.DataFrame:
    """
    For every combination of studies, compute number of shared SNPs,
    and write missing SNPs per study to text files.
    """
    studies = df['study'].unique().to_list()
    for s in studies:
        print(
            f"SNPs that exist in {s} and not in the other studies "
            f"are written to {s}_missing_snps_from_all_others.txt"
        )

    combos = []
    for r in range(1, len(studies) + 1):
        combos += list(combinations(studies, r))

    results = []
    for combo in combos:
        row = {st: 'X' if st in combo else 'O' for st in studies}
        common = set(df.filter(pl.col('study') == combo[0])['SNP'].to_list())
        for st in combo[1:]:
            common &= set(df.filter(pl.col('study') == st)['SNP'].to_list())
        row['number_of_SNPs'] = len(common)
        results.append(row)

        # write missing per study
        for st in studies:
            snps_st = set(df.filter(pl.col('study') == st)['SNP'].to_list())
            others = [o for o in studies if o != st]
            union_others = set().union(*[
                set(df.filter(pl.col('study') == o)['SNP'].to_list()) for o in others
            ])
            missing = union_others - snps_st
            with open(f"{st}_missing_snps_from_all_others.txt", 'w') as f:
                f.write("\n".join(missing))

    return pl.DataFrame(results)


def compute_common_snps(
    df: pl.DataFrame
) -> pl.DataFrame:
    """
    For every combination of studies, compute number of shared SNPs,
    and write missing SNPs per study to text files.
    """
    studies = df['study'].unique().to_list()
    for s in studies:
        print(
            f"SNPs that exist in {s} and not in the other studies "
            f"are written to {s}_missing_snps_from_all_others.txt"
        )

    combos = []
    for r in range(1, len(studies) + 1):
        combos += list(combinations(studies, r))

    results = []
    for combo in combos:
        row = {st: 'X' if st in combo else 'O' for st in studies}
        common = set(df.filter(pl.col('study') == combo[0])['SNP'].to_list())
        for st in combo[1:]:
            common &= set(df.filter(pl.col('study') == st)['SNP'].to_list())
        row['number_of_SNPs'] = len(common)
        results.append(row)

        # write missing per study
        for st in studies:
            snps_st = set(df.filter(pl.col('study') == st)['SNP'].to_list())
            others = [o for o in studies if o != st]
            union_others = set().union(*[
                set(df.filter(pl.col('study') == o)['SNP'].to_list()) for o in others
            ])
            missing = union_others - snps_st
            with open(f"{st}_missing_snps_from_all_others.txt", 'w') as f:
                f.write("\n".join(missing))

    return pl.DataFrame(results)


def main():

    print("-----------------------------------------------------")
    print("PYRAMA - Quality Control")
    print("-----------------------------------------------------")

    parser = argparse.ArgumentParser(
        description="Combine multiple GWAS files, perform QC, and allele harmonization."
    )
    parser.add_argument(
        '--input_files', nargs='+',
        help='Input tab-delimited GWAS files with SNP, CHR, BP, A1, A2, BETA, SE.'
    )
    parser.add_argument(
        '--output', default='merged_study.txt',
        help='Output file for merged and harmonized GWAS data.'
    )
    parser.add_argument(
        '--skip_harm', action='store_true',
        help='Skip allele harmonization and A1/A2 filtering.'
    )
    args = parser.parse_args()

    snp_counts: dict[str, int] = {}
    dfs: list[pl.DataFrame] = []

    for fp in args.input_files:
        print(f"Processing file: {fp}")
        df = read_gwas_file(fp)
        if df is None:
            print(f"Skipping {fp} due to errors.")
            continue
        study = os.path.basename(fp)
        df = df.with_columns(pl.lit(study).alias('study'))
        snp_counts[study] = df['SNP'].n_unique()
        dfs.append(df)

    if not dfs:
        print("No valid GWAS files provided. Exiting.")
        return

    merged = pl.concat(dfs)
    print(f"Total merged records before QC: {merged.height}")

    merged, init_count, final_count, prob_qc = quality_control(
        merged, skip_harm=args.skip_harm
    )
    merged = merged.sort(['SNP', 'CHR'])

    harm_count = 0
    prob_allele = pl.DataFrame({'study': [], 'SNP': [], 'reason': []})

    if not args.skip_harm and {'A1', 'A2'}.issubset(merged.columns):
        merged, prob_allele = filter_consistent_snps(merged)
        merged, harm_count = harmonize_variants(merged)
        print(f"Total allele harmonizations performed: {harm_count}")
    elif args.skip_harm:
        print("Harmonization skipped as requested.")
    else:
        print("Allele columns not found. Harmonization not performed.")

    print(f"Total merged records after QC: {merged.height}")

    # Write per-study harmonized files
    for study in merged['study'].unique().to_list():
        df_st = merged.filter(pl.col('study') == study).drop('study')
        out_file = f"{study.split('.')[0]}_harmonized.txt"
        df_st.write_csv(out_file, separator='\t')
        print(f"Written {df_st.height} rows to {out_file}")

    # Write merged output
    merged.write_csv(args.output, separator='\t')
    print(f"Merged and harmonized data written to {args.output}")

    # 1. Convert each to pandas
    prob_qc_pd     = prob_qc.to_pandas()
    prob_allele_pd = prob_allele.to_pandas()

    # 2. Concatenate in pandas
    all_prob_pd = pd.concat([prob_qc_pd, prob_allele_pd], ignore_index=True)

    # 3. Write out if non‐empty
    if not all_prob_pd.empty:
        all_prob_pd.to_csv('prob_lines.txt', sep='\t', index=False)
        print("Problematic lines written to prob_lines.txt")

    # Compute & write common SNP combinations
    common_df = compute_common_snps(merged)
    common_df.write_csv('common_SNPs_combinations.txt', separator='\t')
    print("Common SNPs combinations written to common_SNPs_combinations.txt")

    # Generate QC report
    report_lines = [
        "Quality Control Report",
        "-----------------------",
        "",
        "Unique SNP counts per study:"
    ]
    for fn, cnt in snp_counts.items():
        report_lines.append(f"{fn}: {cnt}")
    report_lines.append(f"\nTotal allele harmonizations performed: {harm_count}")
    if init_count != final_count:
        report_lines.append(
            f"\nDropped Lines: {init_count - final_count}. "
            "See prob_lines.txt for details.\n"
        )

    with open('quality_control_report.txt', 'w') as rf:
        rf.write("\n".join(report_lines))
    print("Quality control report written to quality_control_report.txt")


if __name__ == '__main__':
    main()
