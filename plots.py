# Imports

import argparse
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from qqman import qqman


def load_data(path: str) -> pl.DataFrame:
    """Load summary stats with Polars and enforce types."""
    # Read only the needed columns, cast types
    df = (
        pl.read_csv(path, separator='\t', columns=['SNP', 'CHR', 'BP', 'P'])
          .with_columns([
              pl.col('CHR').cast(pl.Int32),
              pl.col('BP').cast(pl.Int32),
              pl.col('P').cast(pl.Float64)
          ])
    )
    return df.to_pandas()


def make_manhattan(df , output: str, dpi: int = 600):
    figure, axes = plt.subplots(figsize=(8, 8))

    qqman.manhattan(df, ax=axes, title="Manhattan plot")
    figure.tight_layout()
    plt.savefig(output+'manhattan_plot', format="png",dpi=dpi)
    plt.clf()
    plt.close()


def make_qq(df , output: str, dpi: int = 600):

    figure, axes =  plt.subplots(figsize=(8, 8))


    qqman.qqplot(df, ax=axes, title="QQ plot" )
    figure.tight_layout()
    plt.savefig(output+'qq_plot', format="png", dpi=dpi)
    plt.clf()
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Generate Manhattan and Q-Q plots (with Polars)"
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Path to input TSV file with SNP,CHR,BP,P'
    )
    parser.add_argument(
        '-o', '--output', default='gwas_plot',
        help='Prefix for output files (default: gwas_plot)'
    )
    parser.add_argument(
        '--dpi', type=int, default=600,
        help='Resolution for figures in dpi (default: 600)'
    )
    args = parser.parse_args()

    # Load and plot
    df = load_data(args.input)
    make_manhattan(df, args.output, dpi=args.dpi)
    make_qq(df, args.output, dpi=args.dpi)
    print(f"Plots saved as {args.output}_manhattan.png and {args.output}_qq.png at {args.dpi} dpi.")


if __name__ == '__main__':
    main()
