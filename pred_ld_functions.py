import dask.dataframe as dd
import numpy as np
import pandas as pd
import gc
import polars as pl 
pd.set_option('display.width', None)


# harmonise funcitons
def get_complement(allele):
    complement_map = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    return complement_map.get(allele, None)


def harmonise_data_TOP_LD(joined_data_ref):
    joined_data_ref["REF_complement"] = joined_data_ref["REF1"].apply(get_complement)
    joined_data_ref["ALT_complement"] = joined_data_ref["ALT1"].apply(get_complement)

    # Modify BETA values
    joined_data_ref["beta_modified"] = joined_data_ref.apply(
        lambda row: row["BETA"] if row["ALT1"] == row["A1"] and row["REF1"] == row["A2"] else
        -row["BETA"] if row["ALT1"] == row["A2"] and row["REF1"] == row["A1"] else
        row["BETA"] if row["ALT_complement"] == row["A1"] and row["REF_complement"] == row["A2"] else
        -row["BETA"] if row["ALT_complement"] == row["A2"] and row["REF_complement"] == row["A1"] else
        None,
        axis=1
    )

    # Remove rows where beta_modified is NULL
    joined_data_ref = joined_data_ref.dropna(subset=["beta_modified"])

    # # Save the full joined data to a file
    # joined_data_ref.to_csv(joined_output_file, sep='\t', index=False, quotechar='"')

    # Rename the columns without subsetting
    # joined_data_ref.columns = joined_data_ref.columns.str.replace("ALT1", "A1").str.replace("REF1", "A2").str.replace("beta_modified", "BETA")

    # Ensure correct data types
    joined_data_ref["SNP"] = joined_data_ref["SNP"].astype(str)
    joined_data_ref["CHR"] = joined_data_ref["CHR"].astype(int)
    #joined_data_ref["BP"] = joined_data_ref["BP"].astype(int)
    joined_data_ref["A1"] = joined_data_ref["ALT1"].astype(str)
    joined_data_ref["A2"] = joined_data_ref["REF1"].astype(str)
    joined_data_ref["BETA"] = joined_data_ref["beta_modified"].astype(float)
    joined_data_ref["SE"] = joined_data_ref["SE"].astype(float)

    return joined_data_ref


def harmonise_data_Pheno_Scanner(joined_data_ref):
    joined_data_ref["REF_complement"] = joined_data_ref["REF2"].apply(get_complement)
    joined_data_ref["ALT_complement"] = joined_data_ref["ALT2"].apply(get_complement)

    # Modify BETA values
    joined_data_ref["beta_modified"] = joined_data_ref.apply(
        lambda row: row["BETA"] if row["ALT2"] == row["A1"] and row["REF2"] == row["A2"] else
        -row["BETA"] if row["ALT2"] == row["A2"] and row["REF2"] == row["A1"] else
        row["BETA"] if row["ALT_complement"] == row["A1"] and row["REF_complement"] == row["A2"] else
        -row["BETA"] if row["ALT_complement"] == row["A2"] and row["REF_complement"] == row["A1"] else
        None,
        axis=1
    )

    # Remove rows where beta_modified is NULL
    joined_data_ref = joined_data_ref.dropna(subset=["beta_modified"])

    # # Save the full joined data to a file
    # joined_data_ref.to_csv(joined_output_file, sep='\t', index=False, quotechar='"')

    # Rename the columns without subsetting
    # joined_data_ref.columns = joined_data_ref.columns.str.replace("ALT1", "A1").str.replace("REF1", "A2").str.replace("beta_modified", "BETA")

    # Ensure correct data types
    joined_data_ref["SNP"] = joined_data_ref["SNP"].astype(str)
    joined_data_ref["CHR"] = joined_data_ref["CHR"].astype(int)
    #joined_data_ref["BP"] = joined_data_ref["BP"].astype(int)
    joined_data_ref["A1"] = joined_data_ref["ALT2"].astype(str)
    joined_data_ref["A2"] = joined_data_ref["REF2"].astype(str)
    joined_data_ref["BETA"] = joined_data_ref["beta_modified"].astype(float)
    joined_data_ref["SE"] = joined_data_ref["SE"].astype(float)

    return joined_data_ref


def harmonise_data_Hap_Map(joined_data_ref):
    joined_data_ref["REF_complement"] = joined_data_ref["REF1"].apply(get_complement)
    joined_data_ref["ALT_complement"] = joined_data_ref["ALT1"].apply(get_complement)

    # Modify BETA values
    joined_data_ref["beta_modified"] = joined_data_ref.apply(
        lambda row: row["BETA"] if row["ALT1"] == row["A1"] and row["REF1"] == row["A2"] else
        -row["BETA"] if row["ALT1"] == row["A2"] and row["REF1"] == row["A1"] else
        row["BETA"] if row["ALT_complement"] == row["A1"] and row["REF_complement"] == row["A2"] else
        -row["BETA"] if row["ALT_complement"] == row["A2"] and row["REF_complement"] == row["A1"] else
        None,
        axis=1
    )

    # Remove rows where beta_modified is NULL
    joined_data_ref = joined_data_ref.dropna(subset=["beta_modified"])

    # # Save the full joined data to a file
    # joined_data_ref.to_csv(joined_output_file, sep='\t', index=False, quotechar='"')

    # Rename the columns without subsetting
    # joined_data_ref.columns = joined_data_ref.columns.str.replace("ALT1", "A1").str.replace("REF1", "A2").str.replace("beta_modified", "BETA")

    # Ensure correct data types
    joined_data_ref["SNP"] = joined_data_ref["SNP"].astype(str)
    joined_data_ref["CHR"] = joined_data_ref["CHR"].astype(int)
    #joined_data_ref["BP"] = joined_data_ref["BP"].astype(int)
    joined_data_ref["A1"] = joined_data_ref["ALT1"].astype(str)
    joined_data_ref["A2"] = joined_data_ref["REF1"].astype(str)
    joined_data_ref["BETA"] = joined_data_ref["beta_modified"].astype(float)
    joined_data_ref["SE"] = joined_data_ref["SE"].astype(float)

    return joined_data_ref


def Hap_Map_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")

    # if population not in ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK', 'CHD', 'GIH', "TSI", 'MEX', "ASW"]:
        # print("This population is not available in HapMap files. Please select a different population...")
        # exit()
    # maf_file = f'ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    # ld_file = f'ref/Hap_Map/ld_chr{chrom}_{population}.parquet'
    # maf_df = pd.read_parquet(maf_file  )
    # print(maf_df)
    # print(maf_df.columns)

    # # Calculating the Minor Allele Frequency (MAF)
    # # maf_df['MAF'] = maf_df[['refallele_freq', 'otherallele_freq']].min(axis=1)

    # maf_df['MAF'] = maf_df['otherallele_freq']

    # # Renaming columns in the DataFrame
    # maf_df = maf_df.rename(columns={'refallele': 'REF', 'otherallele': 'ALT'})

    # # Filter MAF DataFrame using Dask
    # # maf_df = dd.read_csv(maf_file, sep='\s+',blocksize=None, usecols=['rs#', 'chrom', 'BP', 'otherallele_freq'],  )

    # # maf_df['het-freq'] = maf_df['het-freq'].astype(float)
    # maf_df = maf_df[maf_df['MAF'] >= float(maf_threshold)]

    # # Process LD DataFrame using Dask
    # ld_df = pd.read_parquet(ld_file,    header=None)

    # ld_df = ld_df[(ld_df[3] != ld_df[4]) & (ld_df[6] >= R2_threshold)]

    # maf_df = maf_df.rename(columns={'rs#': 'rsID'})

    # # Define the new column names
    # new_column_names = {
        # 0: 'pos1',
        # 1: 'pos2',
        # 2: 'pop',
        # 3: 'rsID1',
        # 4: 'rsID2',
        # 5: 'Dprime',
        # 6: 'R2'
    # }

    # # Rename the columns
    # ld_df = ld_df.rename(columns=new_column_names)

    # merged_df = dd.merge(ld_df, maf_df, left_on='rsID1', right_on='rsID')
    # merged_df = dd.merge(merged_df, maf_df, left_on='rsID2', right_on='rsID')

    # merged_df = merged_df[
        # ['pos1', 'pos2', 'rsID1', 'rsID2', 'MAF_x', "MAF_y", 'REF_x', 'REF_y', 'ALT_x', 'ALT_y', "R2", "Dprime"]]
    # merged_df = merged_df.rename(
        # columns={'MAF_x': 'MAF1', 'MAF_y': 'MAF2', 'REF_x': 'REF1', 'REF_y': 'REF2', 'ALT_x': 'ALT1', 'ALT_y': 'ALT2'})

    # if imp_snp_list:
        # final_result = merged_df[merged_df['rsID1'].isin(rs_list) & merged_df['rsID2'].isin(imp_snp_list)]

    # else:
        # final_result = merged_df[merged_df['rsID1'].isin(rs_list)]
    # # Clean up
    # del ld_df, maf_df,merged_df
    # gc.collect()
    # #final_result = final_result.compute()  # Important: This triggers the actual computation
    # if final_result.empty:
        # print("No SNPs found")
        # #
        # #

    # final_result.reset_index(inplace=True, drop=True)
    # #final_result.to_csv('LD_info_Hap_Map_chr' + str(chrom) + '.txt', sep="\t", index=False)

    # # final_result.rename(columns={"MAF1": "ALT_AF1", "MAF2": "ALT_AF2"}).to_csv(
        # # 'LD_info_Hap_Map_chr' + str(chrom) + '.txt', sep="\t", index=False
    # # )
    # final_result.rename(columns={"MAF1": "ALT_AF1", "MAF2": "ALT_AF2"})
    # return final_result
    # 1) Validate population
    valid_pops = ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK',
                  'CHD', 'GIH', 'TSI', 'MEX', 'ASW']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    ld_file  = f'ref/Hap_Map/ld_chr{chrom}_{population}.parquet'

    # 3) Read and preprocess MAF table
    maf_df = (
        pl.read_parquet(maf_file)
          # compute MAF as 'otherallele_freq'
          .with_columns(pl.col("otherallele_freq").alias("MAF"))
          # keep only the needed columns
          .select(["rs#", "refallele", "otherallele", "MAF"])
          # rename to match downstream joins
          .rename({
              "rs#": "rsID",
              "refallele": "REF",
              "otherallele": "ALT"
          })
          # filter on MAF threshold
          .filter(pl.col("MAF") >= float(maf_threshold))
    )

    # 4) Read and preprocess LD table
    #    Parquet likely comes with integer column names (0,1,2,…). We filter & rename by position.
    ld_raw = pl.read_parquet(ld_file)
    cols = ld_raw.columns  # e.g. ['0','1','2',…] or positional ints

    ld_df = (
        ld_raw
        # keep only pairs where alleles differ and R2 >= threshold
        .filter(
            (pl.col(cols[3]) != pl.col(cols[4])) &
            (pl.col(cols[6]) >= float(R2_threshold))
        )
        # rename the first 7 columns to names
        .rename({
            cols[0]: "pos1",
            cols[1]: "pos2",
            cols[2]: "pop",
            cols[3]: "rsID1",
            cols[4]: "rsID2",
            cols[5]: "Dprime",
            cols[6]: "R2"
        })
        # drop any extra columns beyond the 7 we care about
        .select(["pos1", "pos2", "pop", "rsID1", "rsID2", "Dprime", "R2"])
    )

    # 5) Join LD ↔ MAF on rsID1 and rsID2
    #    First join on rsID1 (adds REF, ALT, MAF from maf_df)
    tmp = ld_df.join(maf_df, left_on="rsID1", right_on="rsID", how="inner")
    #    Then join that result on rsID2 (suffix="_2" for the second maf columns)
    merged = tmp.join(
        maf_df,
        left_on="rsID2",
        right_on="rsID",
        how="inner",
        suffix="_2"
    )

    # 6) Select & reorder columns, and rename duplicates
    result = merged.select([
        pl.col("pos1"),
        pl.col("pos2"),
        pl.col("rsID1"),
        pl.col("rsID2"),
        pl.col("MAF").alias("MAF1"),
        pl.col("MAF_2").alias("MAF2"),
        pl.col("REF").alias("REF1"),
        pl.col("REF_2").alias("REF2"),
        pl.col("ALT").alias("ALT1"),
        pl.col("ALT_2").alias("ALT2"),
        pl.col("R2"),
        pl.col("Dprime")
    ])

    # 7) Apply the user’s SNP‐list filter
    if imp_snp_list:
        result = result.filter(
            pl.col("rsID1").is_in(rs_list) &
            pl.col("rsID2").is_in(imp_snp_list)
        )
    else:
        result = result.filter(pl.col("rsID1").is_in(rs_list))

    # # 8) Rename MAF columns to ALT_AF and convert to Pandas
    # final_pl = result.rename({
        # "MAF1": "ALT_AF1",
        # "MAF2": "ALT_AF2"
    # })

    # Convert to Pandas and reset index
    final_pd = result.to_pandas().reset_index(drop=True)

    # 9) Clean up
    del ld_raw, ld_df, maf_df, tmp, merged, result 
    gc.collect()

    if final_pd.empty:
        print("No SNPs found")

    return final_pd

def Hap_Map_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = Hap_Map_LD_info_dask(list(study_df['SNP']), chromosome, population, maf_input, r2threshold,
                                      imp_snp_list)

    outputData = pd.merge(outputData, study_df, left_on='rsID1', right_on='SNP', how='left')
    outputData['CHR'] = chromosome
    # print(outputData.head())

    # Harmonise data with Hap Map
    print("Harmonise data with Hap Map...")
    outputData = harmonise_data_Hap_Map(outputData)

    outputData = outputData.drop(['SNP'], axis=1)

    outputData['imputed'] = 0

    outputData = outputData.groupby('rsID2').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)

    out_df = pd.DataFrame({
        'SNP': outputData['rsID2'],
        'CHR': outputData['CHR'],
        'BP': outputData['pos2'],
        'A1': outputData['ALT2'],
        'A2': outputData['REF2'],
        'BETA': outputData['BETA'],
        'SE': outputData['SE'],

    })

    # OLD
    # pa = 1 - outputData['MAF1'].astype(float)
    # pb = 1 - outputData['MAF2'].astype(float)
    # pA = outputData['MAF1'].astype(float)
    # pB = outputData['MAF2'].astype(float)
    # pT = 1 - outputData['MAF1'].astype(float)
    # pM = 1 - outputData['MAF2'].astype(float)
    #
    # Update the values
    outputData['MAF1'] = outputData['MAF1'].apply(lambda x: 0.4999 if x == 0.5 else x)
    outputData['MAF2'] = outputData['MAF2'].apply(lambda x: 0.4999 if x == 0.5 else x)

    pa = 1 - outputData['MAF1'].astype(float)
    pb = 1 - outputData['MAF2'].astype(float)
    pA = outputData['MAF1'].astype(float)
    pB = outputData['MAF2'].astype(float)
    pT = outputData['MAF1'].astype(float)
    pM = outputData['MAF2'].astype(float)

    r2_value = outputData['R2']
    Dprime = outputData['Dprime']
    D = np.sqrt(r2_value * (pA * pB * pa * pb))

    Dmax_neg_D = np.minimum(pA * pB, (1 - pA) * (1 - pB))
    Dmax_pos_D = np.minimum(pA * (1 - pB), pB * (1 - pA))

    # Calculate the expressions
    expr_neg = abs(Dprime - abs(D / Dmax_neg_D))
    expr_pos = abs(Dprime - abs(D / Dmax_pos_D))

    # Find the minimum expression for each element
    min_expr = pd.DataFrame({'neg': expr_neg, 'BP': expr_pos}).min(axis=1)

    # Find minimum discrepancy and update D
    condition = expr_neg <= expr_pos
    D[condition] = -D[condition]

    outputData['BETA'] = np.where(D < 0, -outputData['BETA'], outputData['BETA'])
    # outputData['BETA'] = np.where (D<0, -outputData['BETA'],outputData['BETA'])
    #
    D = np.sqrt(r2_value * (pA * pB * pa * pb))

    OR_t = np.exp(outputData['BETA'])
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))

    # Sympy

    # f_deriv = (-D * (-D + pT * (1 - pM)) * (OR_t - 1) * OR_t / (
    #             pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1) ** 2) + D * OR_t / (
    #                        pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1))) / (
    #                       D * (OR_t - 1) / (pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1)) + 1)

    beta_t = outputData['BETA']
    numerator = (
            (-D * np.exp(beta_t) * (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM))) /
            ((1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) ** 2 * pM)
            + (D * np.exp(beta_t)) /
            ((1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) * pM)
    )
    denominator = 1 + (D * (-1 + np.exp(beta_t))) / (
            (1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) * pM
    )

    f_deriv = numerator / denominator

    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['BETA'] = np.log(OR_m)
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['BETA'] / out_df['SE']
    out_df['imputed'] = 1
    out_df['R2'] = r2_value
    print(f"Imputed : {sum(out_df['imputed'])} SNPs in chromosome " + str(chromosome))
    # print(out_df.head())
    return out_df


def pheno_Scanner_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
     # 1) Enforce minimum R2 threshold
    if R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 ≥ 0.8. Setting R2_threshold = 0.8")
        R2_threshold = 0.8

    print("Building lazy plan for Pheno Scanner files...")

    # 2) Paths & pop mapping
    maf_file = "ref/Pheno_Scanner/1000G.parquet"
    ld_file  = f"ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.parquet"
    pop_map  = {"EUR":"eur","EAS":"eas","AFR":"afr","AMR":"amr","SAS":"sas"}
    maf_pop  = pop_map.get(population)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # 3) First MAF scan (for ref side → MAF1)
    maf1_lazy = (
        pl.scan_parquet(maf_file)
          .select(["hg19_coordinates","chr","rsid", maf_pop, "a1","a2"])
          .filter(
              (pl.col("chr")==chrom) &
              (pl.col(maf_pop) != "-")
          )
          .with_columns(pl.col(maf_pop).cast(pl.Float64))
          .filter(pl.col(maf_pop) >= maf_threshold)
          .rename({
              "hg19_coordinates":"ref_hg19_coordinates",
              "rsid":"ref_rsid",
              maf_pop:"MAF1",
              "a1":"ALT1",
              "a2":"REF1",
          })
          .select(["ref_hg19_coordinates","ref_rsid","MAF1","ALT1","REF1"])
    )

    # 4) LD scan
    ld_lazy = (
        pl.scan_parquet(ld_file)
          .select(["ref_hg19_coordinates","ref_rsid","rsid","r2","r","dprime"])
          .filter(
              (pl.col("ref_rsid") != pl.col("rsid")) &
              (pl.col("r2") >= R2_threshold)
          )
    )

    # 5) Second MAF scan (for query side → MAF2), keep `rsid` for join!
    maf2_lazy = (
        pl.scan_parquet(maf_file)
          .select(["hg19_coordinates","chr","rsid", maf_pop, "a1","a2"])
          .filter(
              (pl.col("chr")==chrom) &
              (pl.col(maf_pop) != "-")
          )
          .with_columns(pl.col(maf_pop).cast(pl.Float64))
          .filter(pl.col(maf_pop) >= maf_threshold)
          .rename({
              "hg19_coordinates":"pos2(hg19)",
              maf_pop:"MAF2",
              "a1":"ALT2",
              "a2":"REF2",
          })
          # note: we keep `rsid` here so we can join on it
          .select(["pos2(hg19)","rsid","MAF2","ALT2","REF2"])
    )

    # 6) Join the three pieces
    joined = (
        ld_lazy
          .join(maf1_lazy, on="ref_rsid", how="inner")
          .join(maf2_lazy, on="rsid",      how="inner")
    )

    # 7) Rename LD columns & the `rsid` from maf2 → final names
    renamed = joined.rename({
        "ref_rsid":"rsID1",
        "rsid":"rsID2",
        "ref_hg19_coordinates":"pos1(hg19)",
        "r2":"R2",
    })

    # 8) Filter by user‐supplied SNP lists
    filtered = renamed.filter(
        pl.col("rsID2").is_in(rs_list) &
        (pl.col("rsID1").is_in(imp_snp_list) if imp_snp_list else pl.lit(True))
    )

    # 9) Extract numeric coords
    with_pos = filtered.with_columns([
        pl.col("pos1(hg19)").str.extract(r".*:(\d+)$",1).cast(pl.Int64),
        pl.col("pos2(hg19)").str.extract(r".*:(\d+)$",1).cast(pl.Int64),
    ])

    # 10) Final projection/order
    final_lazy = with_pos.select([
        "rsID1","pos1(hg19)","rsID2","dprime",
        "pos2(hg19)","R2","r",
        "MAF1","MAF2","ALT1","REF1","ALT2","REF2",
    ])

    # 11) Execute in streaming, collect to Polars → pandas
    final_pl = final_lazy.collect()
    final_pd = final_pl.to_pandas()

    # 12) Cleanup
    gc.collect()
    if final_pd.empty:
        print("No SNPs found")

    return final_pd

def pheno_Scanner_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = pheno_Scanner_LD_info_dask(list(study_df['SNP']), chromosome, population, maf_input, r2threshold,
                                            imp_snp_list)
    outputData = pd.merge(outputData, study_df, left_on='rsID2', right_on='SNP', how='left')

    outputData['CHR'] = chromosome
    ## Harmonise with Pheno Scanner

    print("Harmonise with Pheno Scanner...")
    outputData = harmonise_data_Pheno_Scanner(outputData)

    ###
    outputData = outputData.drop(['SNP'], axis=1)

    outputData = outputData.rename(
        columns={"rsID2": "rsID1", "rsID1": "rsID2", "MAF1": "MAF2", "MAF2": "MAF1", "pos1(hg19)": "pos2(hg19)",
                 "pos2(hg19)": "pos1(hg19)"})

    outputData['imputed'] = 0
    # print(outputData.head())
    outputData = outputData.groupby('rsID2').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)

    out_df = pd.DataFrame({
        'SNP': outputData['rsID2'],
        'CHR': outputData['CHR'],
        'BP': outputData['pos2(hg19)'],
        'A1': outputData['ALT1'],
        'A2': outputData['REF1'],
        'BETA': outputData['BETA'],
        'SE': outputData['SE'],

    })

    r2_value = outputData['R2']
    r_value = outputData['r']
    # outputData['MAF2']  =  np.where(r_value<0, 1-outputData['MAF2'] ,outputData['MAF2'])
    # outputData['MAF2']  =  np.where(r_value<0, 1-outputData['MAF1'] ,outputData['MAF1'])

    pa = 1 - outputData['MAF1'].astype(float)
    pb = 1 - outputData['MAF2'].astype(float)
    pA = outputData['MAF1'].astype(float)
    pB = outputData['MAF2'].astype(float)
    pT = outputData['MAF1'].astype(float)
    pM = outputData['MAF2'].astype(float)

    Dprime = outputData['dprime']
    D = r_value * np.sqrt(pA * pB * pa * pb)

    outputData['BETA'] = np.where(r_value < 0, -outputData['BETA'], outputData['BETA'])

    #
    #     Dmax_neg_D = np.minimum(pA * pB, (1 - pA) * (1 - pB))
    #     Dmax_pos_D = np.minimum(pA * (1 - pB), (1 - pA) * pB)
    #
    #     # Calculate the expressions
    #     expr_neg = abs(Dprime - abs(D / Dmax_neg_D))
    #     expr_pos = abs(Dprime - abs(D / Dmax_pos_D))
    #
    #     # Find the minimum expression for each element
    #     min_expr = pd.DataFrame({'neg': expr_neg, 'BP': expr_pos}).min(axis=1)
    #
    #     # Check which elements meet the condition
    #     condition = min_expr == expr_neg
    #
    #     # Update D where the condition is true
    #     D[condition] = -D[condition]

    OR_t = np.exp(outputData['BETA'])
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))

    # OLD
    # f_deriv = (-D * (-D + pT * (1 - pM)) * (OR_t - 1) * OR_t / (
    #             pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1) ** 2) + D * OR_t / (
    #                        pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1))) / (
    #                       D * (OR_t - 1) / (pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1)) + 1)

    beta_t = outputData['BETA']
    numerator = (
            (-D * np.exp(beta_t) * (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM))) /
            ((1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) ** 2 * pM)
            + (D * np.exp(beta_t)) /
            ((1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) * pM)
    )
    denominator = 1 + (D * (-1 + np.exp(beta_t))) / (
            (1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) * pM
    )

    f_deriv = numerator / denominator

    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['BETA'] = np.log(OR_m)
    out_df['BETA'] = np.where(r_value < 0, -out_df['BETA'], out_df['BETA'])
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['BETA'] / out_df['SE']
    out_df['imputed'] = 1
    out_df['R2'] = r2_value

    print(f"Imputed : {sum(out_df['imputed'])} SNPs in chromosome " + str(chromosome))
    # print(out_df.head())
    return out_df


 

def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None):
    print("Loading TOP-LD files...")
    
     # 1) Combine rsIDs up‐front
    if imp_snp_list:
        all_rsids = list(set(rs_list) | set(imp_snp_list))
    else:
        all_rsids = rs_list

    # 2) Lazy‐scan the MAF file, project only needed cols, then filter
    maf_path = (
        f"ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    )
    maf_lazy = (
        pl.scan_parquet(maf_path)
          .select(["Uniq_ID", "rsID", "MAF", "REF", "ALT"])
          .filter(pl.col("MAF") >= maf_threshold)
          .filter(pl.col("rsID").is_in(all_rsids))
    )

    # 3) Lazy‐scan the LD file, project and filter by R²
    ld_path = (
        f"ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"
    )
    ld_lazy = (
        pl.scan_parquet(ld_path)
          .select(["Uniq_ID_1", "Uniq_ID_2", "R2", "+/-corr", "Dprime"])
          .filter(pl.col("R2") >= R2_threshold)
    )

    # 4) Prepare two MAF views for joining
    maf1 = maf_lazy.rename({
        "Uniq_ID": "Uniq_ID_1", "rsID": "rsID1", "MAF": "MAF1",
        "REF":      "REF1", "ALT": "ALT1"
    })
    maf2 = maf_lazy.rename({
        "Uniq_ID": "Uniq_ID_2", "rsID": "rsID2", "MAF": "MAF2",
        "REF":      "REF2", "ALT": "ALT2"
    })

    # 5) Join LD ↔ MAF1 ↔ MAF2 all lazily
    joined = (
        ld_lazy
          .join(maf1, on="Uniq_ID_1", how="inner")
          .join(maf2, on="Uniq_ID_2", how="inner")
    )

    # 6) If you provided an imp_snp_list, filter rsID roles
    if imp_snp_list:
        joined = joined.filter(
            (pl.col("rsID1").is_in(rs_list)) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )

    # 7) Select + rename final output columns
    final_lazy = joined.select([
        pl.col("Uniq_ID_1").alias("pos1"),
        pl.col("Uniq_ID_2").alias("pos2"),
        "R2", "+/-corr", "Dprime",
        "rsID1", "rsID2",
        "MAF1", "MAF2",
        "REF1", "ALT1", "REF2", "ALT2"
    ])

    # 8) Execute once in streaming mode (__very__ memory‐efficient)
   
    result = final_lazy.collect()

    if result.is_empty():
        print("No SNPs found.")

    # 9) Cleanup
    del maf_lazy, maf1, maf2, ld_lazy, joined, final_lazy
    gc.collect()

    return result.to_pandas()

def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data

    outputData = TOP_LD_info(list(study_df['SNP']), chromosome, population, maf_input, r2threshold, imp_snp_list)
    # Harmonise with TOP-LD
    print("Harmonise with TOP-LD panel...")

    # study_df = harmonise_data_TOP_LD(study_df,outputData)
    ##
    outputData = pd.merge(outputData, study_df, left_on='rsID1', right_on='SNP', how='left')
    outputData = harmonise_data_TOP_LD(outputData)

    outputData = outputData.drop(['SNP', 'BP'], axis=1)

    outputData['imputed'] = 0
    outputData = outputData.groupby('rsID2').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)

    out_df = pd.DataFrame({
        'SNP': outputData['rsID2'],
        'CHR': outputData['CHR'],
        'BP': outputData['pos2'],
        'A1': outputData['ALT2'],
        'A2': outputData['REF2'],
        'BETA': outputData['BETA'],
        'SE': outputData['SE'],

    })
    # # ##OLD###
    # pa = 1 - outputData['MAF1']
    # pb = 1 - outputData['MAF2']
    # pA = outputData['MAF1']
    # pB = outputData['MAF2']
    # pT = 1 - outputData['MAF1']
    # pM = 1 - outputData['MAF2']
    # # ####
    # #
    #

    # ##NEW###
    pa = 1 - outputData['MAF1']
    pb = 1 - outputData['MAF2']
    pA = outputData['MAF1']
    pB = outputData['MAF2']
    pT = outputData['MAF1']
    pM = outputData['MAF2']
    # ####
    #

    Dprime = outputData['Dprime']
    r2_value = outputData['R2']
    sign = outputData['+/-corr']

    r_value = np.sqrt(r2_value)
    r_value = np.where(sign == '-', -r_value, r_value)
    D = np.sqrt(r2_value * (pA * pB * pa * pb))

    # # Change the sign of D if sign is '-'
    # D = np.where(sign == '-', -D, D)
    #
    # outputData['BETA']  =  np.where(r_value<0, -outputData['BETA'] ,outputData['BETA'])

    #
    #
    # Dmax_neg_D = np.minimum(pa*pb, (1-pa)*(1-pb))
    # Dmax_pos_D = np.minimum(pa*(1-pb), (1-pa)*pb)
    # #
    # # # Calculate the expressions
    # expr_neg = abs(Dprime - abs(D / Dmax_neg_D))
    # expr_pos = abs(Dprime - abs(D / Dmax_pos_D))
    # #
    # # ## Find the minimum expression for each element
    # min_expr = pd.DataFrame({'neg': expr_neg, 'BP': expr_pos}).min(axis=1)
    # #
    # # # Check which elements meet the condition
    # condition = min_expr == expr_neg
    # #
    # # #Update D where the condition is true
    # D[condition] = -D[condition]

    outputData['BETA'] = np.where(r_value < 0, -outputData['BETA'], outputData['BETA'])

    OR_t = np.exp(outputData['BETA'])

    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))

    # Sympy OLD
    # f_deriv = (-D * (-D + pT * (1 - pM)) * (OR_t - 1) * OR_t / (
    #             pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1) ** 2) + D * OR_t / (
    #                        pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1))) / (
    #                       D * (OR_t - 1) / (pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1)) + 1)
    #

    beta_t = outputData['BETA']
    numerator = (
            (-D * np.exp(beta_t) * (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM))) /
            ((1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) ** 2 * pM)
            + (D * np.exp(beta_t)) /
            ((1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) * pM)
    )
    denominator = 1 + (D * (-1 + np.exp(beta_t))) / (
            (1 + (-1 + np.exp(beta_t)) * (-D + pT * (1 - pM)) - pM) * pM
    )

    f_deriv = numerator / denominator

    Var_M = f_deriv ** 2 * var_x ** 2
    # Var_M = (1/OR_m) ** 2 * var_x ** 2

    out_df['BETA'] = np.log(OR_m)

    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['BETA'] / out_df['SE']
    out_df['imputed'] = 1
    out_df['R2'] = r2_value

    print(f"Imputed : {sum(out_df['imputed'])} SNPs in chromosome " + str(chromosome))
    
    return out_df


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    final_results_list = []

    study_df = pd.read_csv(file_path, sep="\t")
    chroms = list(set(study_df['CHR']))


    # Take only 22 chromosomes
    chroms= [chrom for chrom in chroms if not str(chrom).isdigit() or int(chrom) <= 22]
    ref_panel = ref_file
    # Check if all required columns are present
    required_columns = ['SNP', 'CHR', 'BP', 'BETA', 'SE']  # The required order of columns

    imputed_columns = [col for col in required_columns if col not in required_columns]
    if imputed_columns:
        raise ValueError(
            f" Warning: Check the column names! The columns must be in the following order: SNP, CHR, BP, BETA, SE")

    # Depending on the reference panel...
    if ref_panel == 'TOP_LD':
        for chrom in chroms:
            final_data = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            data = pd.read_csv(file_path, sep="\t")

            data['z'] = data['BETA'] / data['SE']
            data['imputed'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            #final_data.to_csv(file_path+"_imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            #print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            #print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'SNP' column.
            final_df_init = final_df_init.drop_duplicates(subset="SNP")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            # Step 1: Identify SNPs that exist both as imputed (1) and non-imputed (0)
            duplicate_snps = final_data[final_data['SNP'].duplicated(keep=False)]['SNP'].unique()
            
            # Step 2: Filter the data
            # Keep non-imputed rows for duplicate SNPs and keep all other rows as they are
            final_data = final_data[~final_data['SNP'].isin(duplicate_snps) | (final_data['imputed'] == 0)]
            
            final_data.to_csv(file_path+"_imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Pheno_Scanner':
        for chrom in chroms:
            final_data = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['BETA'] / data['SE']
            data['imputed'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

           # print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            #print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'SNP' column.
            final_df_init = final_df_init.drop_duplicates(subset="SNP")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            # Step 1: Identify SNPs that exist both as imputed (1) and non-imputed (0)
            duplicate_snps = final_data[final_data['SNP'].duplicated(keep=False)]['SNP'].unique()
            
            # Step 2: Filter the data
            # Keep non-imputed rows for duplicate SNPs and keep all other rows as they are
            final_data = final_data[~final_data['SNP'].isin(duplicate_snps) | (final_data['imputed'] == 0)]
            
            final_data.to_csv(file_path+"_imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Hap_Map':
        for chrom in chroms:
            final_data = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['BETA'] / data['SE']
            data['imputed'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            #print(f"Total : {len(final_data)} SNPs")

            #final_data.to_csv(file_path+"_imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

           # print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
           # print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'SNP' column.
            final_df_init = final_df_init.drop_duplicates(subset="SNP")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            
            # Step 1: Identify SNPs that exist both as imputed (1) and non-imputed (0)
            duplicate_snps = final_data[final_data['SNP'].duplicated(keep=False)]['SNP'].unique()
            
            # Step 2: Filter the data
            # Keep non-imputed rows for duplicate SNPs and keep all other rows as they are
            final_data = final_data[~final_data['SNP'].isin(duplicate_snps) | (final_data['imputed'] == 0)]
            
            final_data.to_csv(file_path+"_imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'all_panels':

        print(f"Checking all LD sources")
        for chrom in chroms:
            # For HapMap, we need to take all the panels and merge them...
            if population == 'EUR':
                pop_hm = "CEU"
                final_data_hm = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                # pop_hm = "TSI"
                # final_data_hm_TSI = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                #

                # final_data_hm = pd.concat([final_data_hm_CEU, final_data_hm_TSI])
                # Keep the largest R2 value if a SNP is common in any of the panels
                final_data_hm = final_data_hm.groupby('SNP').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(
                    drop=True)

            if population == 'AFR':
                pop_hm = "YRI"
                final_data_hm_YRI = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                pop_hm = "MKK"
                final_data_hm_MKK = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                pop_hm = "LWK"
                final_data_hm_LWK = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                pop_hm = "ASW"
                final_data_hm_ASW = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                final_data_hm = pd.concat([final_data_hm_YRI, final_data_hm_MKK, final_data_hm_LWK, final_data_hm_ASW])

            if population == 'EAS':
                pop_hm = "CHB"
                final_data_hm_CHB = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)
                pop_hm = "JPT"
                final_data_hm_JPT = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)
                pop_hm = "CHD"
                final_data_hm_CHD = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                final_data_hm = pd.concat([final_data_hm_CHB, final_data_hm_JPT, final_data_hm_CHD])

            if population == 'SAS':
                pop_hm = 'GIH'
                final_data_hm = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

            # Keep the largest R2 value if a SNP is common in any of the panels
            final_data_hm = final_data_hm.groupby('SNP').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
            final_data_hm['source'] = 'HapMap'

            # final_data_hm = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_ps = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_tld = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_ps['source'] = 'Pheno Scanner'
            final_data_tld['source'] = 'TOP-LD'

            # final_data = pd.concat([final_data_hm, final_data_ps, final_data_tld])
            final_data = pd.concat([final_data_hm, final_data_ps, final_data_tld], ignore_index=True)

            # Keep the largest R2 value if a SNP is common in any of the panels
            if imp_snp_list == True:
                final_data = final_data.loc[final_data.groupby('SNP')['R2'].idxmax()]
            else:
                final_data = final_data.groupby('SNP').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
            #print(len(final_data))

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['BETA'] / data['SE']
            data['imputed'] = 0
            data['source'] = 'GWAS'
            final_data = pd.concat([final_data, data], ignore_index=True)

            print(f"Total Imputed SNPs: {len(final_data[final_data['imputed'] == 1])} SNPs")

            #print(f"Total : {len(final_data)} SNPs")

            #final_data.to_csv(file_path+"_imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            #print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            #print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'SNP' column.
            final_df_init = final_df_init.drop_duplicates(subset="SNP")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            
            # Step 1: Identify SNPs that exist both as imputed (1) and non-imputed (0)
            duplicate_snps = final_data[final_data['SNP'].duplicated(keep=False)]['SNP'].unique()
            
            # Step 2: Filter the data
            # Keep non-imputed rows for duplicate SNPs and keep all other rows as they are
            final_data = final_data[~final_data['SNP'].isin(duplicate_snps) | (final_data['imputed'] == 0)]
            
            final_data.to_csv(file_path+"_imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")
