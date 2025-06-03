import pandas as pd
import os
from glob import glob

def combine_trait_files(folder_path):
    """
    Combines multiple trait files in a given folder. Each file should contain columns for 
    'snp', 'chr', 'pos', 'beta', 'SE', 'z', and 'p'. A 'trait' column will be added to 
    the final DataFrame based on the filename of each file.

    Parameters:
    folder_path (str): The path to the folder containing the trait files.

    Returns:
    pd.DataFrame: A combined DataFrame with an additional 'trait' column.
    """
    # Define the columns we want to keep
    columns_of_interest = ['snp', 'chr', 'pos', 'beta', 'SE', 'z', 'p']
    
    # Initialize an empty list to hold DataFrames
    df_list = []
    
    # Get all CSV/TSV files in the folder
    file_paths = glob(os.path.join(folder_path, "*.csv")) + glob(os.path.join(folder_path, "*.tsv"))
    
    for file_path in file_paths:
        # Extract the filename (without extension) as the trait name
        trait_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # Read the file into a DataFrame (automatically handles CSV and TSV)
        df = pd.read_csv(file_path, sep=None, engine='python')
        
        # Filter for columns of interest, handling missing columns
        df = df.loc[:, [col for col in columns_of_interest if col in df.columns]]
        
        # Add the 'trait' column
        df['trait'] = trait_name
        
        # Append the DataFrame to the list
        df_list.append(df)
    
    # Combine all DataFrames into one
    combined_df = pd.concat(df_list, ignore_index=True)
    
    return combined_df

# Example usage:
# folder_path = "/path/to/trait/files"
# combined_data = combine_trait_files(folder_path)
# print(combined_data)
