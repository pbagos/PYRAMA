 #Imports
import argparse
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from gprofiler import GProfiler

def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform g:Profiler enrichment analysis on SNPs"
    )
    parser.add_argument('-i', '--input', required=True,
                        help="Input tab-delimited file with 'SNP' and 'P' columns")
    parser.add_argument('-o', '--output', default='gprofiler_results.csv',
                        help="Output CSV file for enrichment results")
    parser.add_argument('-p', '--plot', default='gprofiler_manhattan.png',
                        help="Output PNG file for Manhattan-style plot")
    parser.add_argument('-s', '--significance', type=float, default=0.05,
                        help="Adjusted p-value significance threshold")
    parser.add_argument('--organism', default='hsapiens',
                        help="Organism name for g:Profiler (e.g. 'hsapiens')")
    parser.add_argument('--sources', nargs='+', default=["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP"],
                        help="Data sources to query (e.g. GO:BP GO:CC KEGG). If omitted, all available sources will be used.")
    parser.add_argument('--no_iea', action='store_true',
                        help="Exclude electronic GO annotations (IEA)")
    return parser.parse_args()

def load_sig_snps(input_file: str, p_threshold: float) -> list:
    df = pl.read_csv(input_file, separator='\t')
    for col in ('SNP', 'P'):
        if col not in df.columns:
            raise ValueError("Input file must contain 'SNP' and 'P' columns")
    sig_df = df.filter(pl.col('P') <= p_threshold)
    snps = sig_df.select('SNP').drop_nulls().to_series().cast(str).to_list()
    if not snps:
        raise ValueError(f"No SNPs found with P <= {p_threshold}")
    return snps

def run_enrichment(snps: list, organism: str, sources: list, exclude_iea: bool) -> pl.DataFrame:
    gp = GProfiler(return_dataframe=True)
    pandas_df = gp.profile(
        organism=organism,
        query=snps,
        sources=sources,
        no_iea=exclude_iea
    )
    if pandas_df is None or pandas_df.empty:
        raise RuntimeError("g:Profiler returned no enriched terms.")
    # Convert any list columns to comma-separated strings
    for col in ('intersection', 'parents'):
        if col in pandas_df.columns:
            pandas_df[col] = pandas_df[col].apply(lambda x: ','.join(x) if isinstance(x, (list, tuple)) else x)
    return pl.from_pandas(pandas_df)

def plot_manhattan(df: pl.DataFrame, significance: float, plot_file: str):

    # convert to pandas and filter by p-value
    df = df.to_pandas()
    df = df[df.p_value < float(significance)].copy()
    df['negative_log10_of_adjusted_p_value'] = -np.log10(df['p_value'])
    #df = df[df['negative_log10_of_adjusted_p_value'] <= 16]

    # assign each source a group_id, sort, and create an index for plotting
    df['group_id'] = df.groupby('source').ngroup()
    df.sort_values(['group_id', 'native'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df['i'] = df.index

    # prepare tick positions and labels from the same grouping
    grouped = (
        df
        .groupby('group_id', as_index=False)
        .agg({
            'i': 'median',       # position
            'source': 'first'    # label
        })
    )

    # make the plot
    plt.figure(dpi=600, figsize=(14, 4))
    plot = sns.scatterplot(
        data=df,
        x='i',
        y='negative_log10_of_adjusted_p_value',
        hue='source',
        palette='bright',
        s=9,
        legend='full'
    )
    plot.axhline(16, linestyle='--', linewidth=1, color='gray')

    # set labels, ticks, and title
    plot.set_xlabel('Source')
    plot.set_ylabel('-log10(Padj)')
    plot.set_title('Enrichment Analysis — Manhattan plot', fontsize=16)

    # now ticks and labels match one‐to‐one
    plot.set_xticks(grouped['i'])
    plot.set_xticklabels(grouped['source'], rotation=25, ha='right')

    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()

def main():
    args = parse_args()
    print(f"Loading SNPs from {args.input} (P <= {args.significance})...")
    snps = load_sig_snps(args.input, args.significance)
    print(f"Running g:Profiler enrichment for {len(snps)} SNPs (organism={args.organism})...")
    res = run_enrichment(
        snps,
        organism=args.organism,
        sources=args.sources,
        exclude_iea=args.no_iea
    )
    sig_res = res.filter(pl.col('p_value') <= args.significance)
    if sig_res.is_empty():
        print("Warning: No enriched terms passed the significance threshold.")
    else:
        print(f"{sig_res.height} enriched terms pass the threshold (p <= {args.significance}).")
    print(f"Saving full results to {args.output}...")
    # Write CSV without nested data issues
    res.write_csv(args.output)
    print(f"Plotting Manhattan plot to {args.plot}...")
    plot_manhattan(sig_res, args.significance, args.plot)
    print("Done.")

if __name__ == '__main__':
    main()
