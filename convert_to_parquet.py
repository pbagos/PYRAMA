import os
import pandas as pd

def detect_separator(filename):
    """Detect separator based on file extension."""
    if filename.endswith(".csv.gz"):
        return ","
    elif filename.endswith(".txt.gz"):
        return "\s+"
    else:
        return None

def convert_gz_to_parquet(root_dir):
    """
    Traverse a directory structure, find all .gz files,
    read them into pandas, and convert to .parquet format.
    """
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.endswith(".csv.gz") or file.endswith(".txt.gz"):
                full_path = os.path.join(dirpath, file)
                separator = detect_separator(file)

                if separator is None:
                    print(f"Skipping file with unknown format: {file}")
                    continue

                try:
                    print(f"Processing: {full_path}")
                    df = pd.read_csv(full_path, sep=separator, compression='gzip')

                    # Construct Parquet filename
                    parquet_filename = file.replace('.csv.gz', '.parquet').replace('.txt.gz', '.parquet')
                    parquet_path = os.path.join(dirpath, parquet_filename)

                    df.to_parquet(parquet_path, engine='pyarrow', index=False)
                    print(f"Saved: {parquet_path}")

                except Exception as e:
                    print(f"Error processing {full_path}: {e}")

# Example usage
if __name__ == "__main__":
    root_directory = "ref/Hap_Map/"
    convert_gz_to_parquet(root_directory)
