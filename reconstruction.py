import sys
import numpy as np
import pandas as pd
import dask.dataframe as dd
from gwas_reconstructions.gwas_reconstruction import gwas_reconstruction


def reconstruction(data_path, output_path ):
    try:
        print("Running Reconstruction:")

        # Read data as Dask DataFrame
        data = dd.read_csv(data_path, sep='\t')
        print("Data loaded successfully.")

        # Calculate the Odds Ratio (OR)
        data['OR'] = np.exp(data['BETA'])


        print("Data computation completed.")

        # Apply GWAS reconstruction
        result = gwas_reconstruction(data, method="direct")
        base_filename = 'output'
        # Save results to a specified file
        #output_path = f"{output_folder}/{base_filename}_reconstruction_results.txt"
        result.to_csv(output_path, sep='\t', index=False)

        print(f"Reconstruction results saved to {output_path}")
        return result

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


if __name__ == "__main__":

    data_path, output_folder = sys.argv[1], sys.argv[2]
    print(f"Data path: {data_path}")
    print(f"Output folder: {output_folder}")


    # Run reconstruction
    reconstruction(data_path, output_folder )
