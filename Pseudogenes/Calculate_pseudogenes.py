import os  # Importing the OS module for interacting with the operating system.
import pandas as pd  # Importing pandas for data manipulation and analysis.
import numpy as np  # Importing numpy for numerical operations.
import gffutils  # Importing gffutils for working with GFF files.
import gzip  # Importing gzip for working with compressed files.

# Function to compute median
def median(lst):
    sorted_lst = sorted(lst)  # Sorting the list.
    length = len(sorted_lst)  # Finding the length of the sorted list.
    if length % 2 == 0:
        return (sorted_lst[length // 2 - 1] + sorted_lst[length // 2]) / 2  # Returning the average of the two middle values if the length is even.
    else:
        return sorted_lst[length // 2]  # Returning the middle value if the length is odd.

# Function to load and filter BLAST results from gzipped TSV files
def load_and_filter_blast(directory):
    transposase_ids = set()  # Initializing a set to store transposase IDs.
    for filename in os.listdir(directory):  # Repeat over files in the given directory.
        if filename.endswith('.tsv.gz'):  # Checking if the file is a gzipped TSV file.
            filepath = os.path.join(directory, filename)  # Constructing the full file path.
            with gzip.open(filepath, 'rt') as f:  # Opening the gzipped file in text mode.
                blast_df = pd.read_csv(f, sep='\t', header=None)  # Reading the file into a pandas DataFrame.
                print(f"Processing {filename} with columns: {blast_df.columns.tolist()}")  # Printing debug information.
                
                description_col = blast_df.columns[-1]  # Getting the last column which contains the description.
                hits = blast_df[blast_df[description_col].str.contains("transposase", case=False, na=False)]  # Filtering rows containing 'transposase'.
                transposase_ids.update(hits[0])  # The first column is 'qseqid' and adding these IDs to the set.
    print(f"Total transposase IDs found: {len(transposase_ids)}")  # Printing the total number of transposase IDs found.
    return transposase_ids  # Returning the set of transposase IDs.

# Function to filter GFF files
def filter_gff(gff_file, exclude_ids):
    with open(gff_file, 'r') as file:  # Opening the GFF file.
        lines = file.readlines()  # Reading all lines from the file.
    filtered_lines = [line for line in lines if not any(ex_id in line for ex_id in exclude_ids)]  # Filtering out lines containing exclude IDs.
    print(f"Filtered {len(lines) - len(filtered_lines)} transposase entries from {gff_file}")  # Printing the number of filtered entries.
    return filtered_lines  # Returning the filtered lines.

# Function to write filtered GFF files
def write_filtered_gff(filtered_data, output_filename):
    with open(output_filename, 'w') as file:  # Opening the output file in write mode.
        file.writelines(filtered_data)  # Writing the filtered data to the file.

# Function to process directories and filter GFF files
def process_directories(blastp_directory, blastx_directory, gff_directory, output_directory):
    transposase_ids = set()  # Initializing a set to store transposase IDs.
    transposase_ids.update(load_and_filter_blast(blastp_directory))  # Adding transposase IDs from BlastP directory.
    transposase_ids.update(load_and_filter_blast(blastx_directory))  # Adding transposase IDs from BlastX directory.
    
    if not os.path.exists(output_directory):  # Checking if the output directory exists.
        os.makedirs(output_directory)  # Creating the output directory if it doesn't exist.

    for filename in os.listdir(gff_directory):  # Iterating over files in the GFF directory.
        if filename.endswith('.gff'):  # Checking if the file is a GFF file.
            gff_file = os.path.join(gff_directory, filename)  # Constructing the full file path.
            filtered_data = filter_gff(gff_file, transposase_ids)  # Filtering the GFF file.
            output_gff_file = os.path.join(output_directory, filename)  # Constructing the output file path.
            write_filtered_gff(filtered_data, output_gff_file)  # Writing the filtered data to the output file.
            print(f'Processed and saved: {output_gff_file}')  # Printing a message indicating the file has been processed and saved.

# Function to categorize genes in GFF files
def categorize_genes(gff_file):
    categories = {  # Initializing categories dictionary.
        'Truncated': [],
        'Run-on': [],
        'Intergenic region': [],
        'Predicted fragmentation': []
    }
    
    with open(gff_file, 'r') as f:  # Opening the GFF file.
        for line in f:  # Iterating over each line in the file.
            if line.startswith('#') or len(line.split('\t')) < 5:
                continue  # Skipping comment lines or lines with fewer than 5 columns.

            start = int(line.split('\t')[3])  # Extracting the start position.
            end = int(line.split('\t')[4])  # Extracting the end position.
            length = end - start + 1  # Calculating the length of the feature.
            
            locus_tag = None
            for attribute in line.split('\t')[8].split(';'):  # Iterating over attributes.
                if attribute.startswith('locus_tag='):
                    locus_tag = attribute.split('=')[1]  # Extracting the locus_tag.
                    break
            
            if not locus_tag:
                continue  # Skipping if locus_tag is not found.

            if 'ORF is' in line:
                orf_percentage = float(line.split('ORF is')[1].split('%')[0].strip())  # Extracting the ORF percentage.
                if orf_percentage < 75:
                    categories['Truncated'].append((locus_tag, length))  # Adding to 'Truncated' category if ORF percentage is less than 75.
                else:
                    categories['Run-on'].append((locus_tag, length))  # Adding to 'Run-on' category if ORF percentage is 75 or more.
            elif 'Intergenic region with' in line:
                categories['Intergenic region'].append((locus_tag, length))  # Adding to 'Intergenic region' category.
            elif 'Predicted fragmentation of a single gene' in line:
                categories['Predicted fragmentation'].append((locus_tag, length))  # Adding to 'Predicted fragmentation' category.
    
    return categories  # Returning the categorized genes.

# Function to categorize genes and compute statistics per category
def categorize_genes_stats(gff_file):
    categories = {  # Initializing categories dictionary with statistics.
        'Truncated': {'count': 0, 'length': 0, 'all_lengths': []},
        'Run-on': {'count': 0, 'length': 0, 'all_lengths': []},
        'Intergenic region': {'count': 0, 'length': 0, 'all_lengths': []},
        'Predicted fragmentation': {'count': 0, 'length': 0, 'all_lengths': []}
    }
    
    with open(gff_file, 'r') as f:  # Opening the GFF file.
        for line in f:  # Iterating over each line in the file.
            if line.startswith('#') or len(line.split('\t')) < 5:
                continue  # Skipping comment lines or lines with fewer than 5 columns.

            start = int(line.split('\t')[3])  # Extracting the start position.
            end = int(line.split('\t')[4])  # Extracting the end position.
            length = end - start + 1  # Calculating the length of the feature.

            if 'ORF is' in line:
                orf_percentage = float(line.split('ORF is')[1].split('%')[0].strip())  # Extracting the ORF percentage.
                if orf_percentage < 75:
                    categories['Truncated']['count'] += 1  # Incrementing the count for 'Truncated' category.
                    categories['Truncated']['length'] += length  # Adding the length to the total length for 'Truncated' category.
                    categories['Truncated']['all_lengths'].append(length)  # Appending the length to the list of lengths for 'Truncated' category.
                else:
                    categories['Run-on']['count'] += 1  # Incrementing the count for 'Run-on' category.
                    categories['Run-on']['length'] += length  # Adding the length to the total length for 'Run-on' category.
                    categories['Run-on']['all_lengths'].append(length)  # Appending the length to the list of lengths for 'Run-on' category.
            elif 'Intergenic region with' in line:
                categories['Intergenic region']['count'] += 1  # Incrementing the count for 'Intergenic region' category.
                categories['Intergenic region']['length'] += length  # Adding the length to the total length for 'Intergenic region' category.
                categories['Intergenic region']['all_lengths'].append(length)  # Appending the length to the list of lengths for 'Intergenic region' category.
            elif 'Predicted fragmentation of a single gene' in line:
                categories['Predicted fragmentation']['count'] += 1  # Incrementing the count for 'Predicted fragmentation' category.
                categories['Predicted fragmentation']['length'] += length  # Adding the length to the total length for 'Predicted fragmentation' category.
                categories['Predicted fragmentation']['all_lengths'].append(length)  # Appending the length to the list of lengths for 'Predicted fragmentation' category.
    
    return categories  # Returning the categorized genes with statistics.

# Function to process GFF files and compute statistics
def process_gff(file_path):
    db = gffutils.create_db(file_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)  # Creating an in-memory database from the GFF file.
    all_features = list(db.all_features())  # Retrieving all features from the database.
    
    feature_count = len(all_features)  # Counting the number of features.
    total_length = sum(len(feature) for feature in all_features)  # Calculating the total length of all features.
    lengths = [len(feature) for feature in all_features]  # Creating a list of lengths of all features.
    median_length = np.median(lengths)  # Calculating the median length of all features.
    
    return feature_count, total_length, median_length  # Returning the feature count, total length, and median length.

# Function to filter out rows with a specific prefix in GFF_File
def filter_out_prefix(input_file, prefix, output_file):
    data = pd.read_csv(input_file, sep="\t")
    filtered_data = data[~data['GFF_File'].str.startswith(prefix)]
    filtered_data.to_csv(output_file, sep="\t", index=False)
    return output_file

# Function to process all GFF files and combine all tasks
def process_all(directory_path, output_combined_file, output_stats_file, output_category_stats_file):
    blastp_directory = os.path.join(directory_path, 'BlastP_Pseudofinder')  # Constructing the path for the BlastP directory.
    blastx_directory = os.path.join(directory_path, 'BlastX_Pseudofinder')  # Constructing the path for the BlastX directory.
    gff_directory = os.path.join(directory_path, 'gff')  # Constructing the path for the GFF directory.
    output_directory = os.path.join(directory_path, 'filtered_gff')  # Constructing the path for the output directory.
    
    # Filter GFF files based on BLAST results
    process_directories(blastp_directory, blastx_directory, gff_directory, output_directory)  # Filtering GFF files based on BLAST results.
    
    # Write the header for the combined file
    with open(output_combined_file, 'w') as out:
        out.write("GFF_File\tCategory\tGene Identifier\tLength\n")  # Writing the header for the combined file.

    # Write the header for the statistics file
    with open(output_stats_file, 'w') as out:
        out.write("GFF_File\tFeature count\tTotal length\tMedian length\n")  # Writing the header for the statistics file.

    # Write the header for the category statistics file
    with open(output_category_stats_file, 'w') as out:
        out.write("GFF_File\tCategory\tCount\tOverall length\tMedian length\n")  # Writing the header for the category statistics file.
    
    for filename in os.listdir(output_directory):  # Iterating over files in the output directory.
        if filename.endswith(".gff"):  # Checking if the file is a GFF file.
            gff_path = os.path.join(output_directory, filename)  # Constructing the full file path.
            
            # Categorize genes and append to the combined file
            categories = categorize_genes(gff_path)  # Categorizing genes.
            with open(output_combined_file, 'a') as out:
                for category, data in categories.items():  # Iterating over categories.
                    for gene_data in data:
                        out.write(f"{filename}\t{category}\t{gene_data[0]}\t{gene_data[1]}\n")  # Writing categorized genes to the combined file.
            
            # Compute statistics and append to the statistics file
            count, total_length, median_length = process_gff(gff_path)  # Computing statistics.
            with open(output_stats_file, 'a') as out:
                out.write(f"{filename}\t{count}\t{total_length}\t{median_length}\n")  # Writing statistics to the statistics file.

            # Compute category statistics and append to the category statistics file
            category_stats = categorize_genes_stats(gff_path)  # Computing category statistics.
            with open(output_category_stats_file, 'a') as out:
                for category, data in category_stats.items():  # Iterating over categories.
                    med_length = median(data['all_lengths']) if data['all_lengths'] else 0  # Calculating the median length for the category.
                    out.write(f"{filename}\t{category}\t{data['count']}\t{data['length']}\t{med_length}\n")  # Writing category statistics to the category statistics file.
            
            print(f'{filename}: Feature count: {count}, Total length: {total_length} bp, Median length: {median_length} bp')  # Printing a summary of the statistics.

    # Filter out rows with the prefix 'filtered_0_genes_removed'
    filter_out_prefix(output_category_stats_file, 'filtered_0_genes_removed', output_category_stats_file)
    print(f"Filtered out rows with prefix 'filtered_0_genes_removed' in {output_category_stats_file}")

# Define the main directory path containing all subdirectories
directory_path = "/Users/vesna/Desktop/Marine_Microbial_Symbiosis/Richelia/Pseudogenes"  # Defining the main directory path.
output_combined_file = os.path.join(directory_path, "combined_results.tsv")  # Defining the path for the combined results file.
output_stats_file = os.path.join(directory_path, "feature_statistics.tsv")  # Defining the path for the feature statistics file.
output_category_stats_file = os.path.join(directory_path, "category_statistics.tsv")  # Defining the path for the category statistics file.

# Run the full processing
process_all(directory_path, output_combined_file, output_stats_file, output_category_stats_file)

