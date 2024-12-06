import os
import pandas as pd

# List of all gc_content.tsv files
files = [f for f in os.listdir() if f.endswith('_gc_content.tsv')]

# List to hold dataframes
dfs = []

for file in files:
    # Extract sample name from filename
    sample = os.path.splitext(file)[0].replace('_gc_content', '')
    
    # Read the file into a dataframe
    df = pd.read_csv(file, sep='\t')
    
    # Add a column for the sample name
    df['Sample'] = sample
    
    # Append to the list of dataframes
    dfs.append(df)

# Concatenate all dataframes
combined_df = pd.concat(dfs, ignore_index=True)

# Reorder columns
combined_df = combined_df[['Sample', 'Category', 'Seq_ID', 'GC_Content']]

# Save the combined dataframe to a new TSV file
combined_df.to_csv('combined_gc_content.tsv', sep='\t', index=False)
