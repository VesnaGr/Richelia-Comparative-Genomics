import os  # Used for file and directory operations
import pandas as pd  # Used to handle and save the analysis results
from statistics import median  # Used to calculate the median of a list of numbers
import gffutils  # Used for GFF file processing

def calculate_total_assembly_size(fasta_file):
    total_size = 0  # Initialize total size to 0
    with open(fasta_file, 'r') as file:  # Open the FASTA file for reading
        for line in file:  # Iterate over each line in the file
            if not line.startswith('>'):  # Ignore header lines that start with '>'
                total_size += len(line.strip())  # Add the length of the sequence line to total size
    return total_size  # Return the total size of the assembly

def process_gff(file_path):
    # Create an in-memory database from the GFF file
    db = gffutils.create_db(file_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    
    cds_sizes = []  # List to store sizes of CDS features
    trna_sizes = []  # List to store sizes of tRNA features
    transposase_sizes = []  # List to store sizes of transposase features
    intergenic_sizes = []  # List to store sizes of intergenic regions
    intergenic_start = None  # Variable to keep track of the end position of the last feature

    for feature in db.all_features(order_by=('seqid', 'start', 'end')):  # Iterate over all features in the GFF file
        feature_type = feature.featuretype  # Get the feature type (e.g., CDS, tRNA, etc.)
        start = feature.start  # Get the start position of the feature
        end = feature.end  # Get the end position of the feature
        size = end - start + 1  # Calculate the size of the feature

        if feature_type == "CDS":  # If the feature is a CDS
            cds_sizes.append(size)  # Append the size to cds_sizes list
            if 'transposase' in feature.attributes.get('product', [''])[0]:  # If the product attribute contains 'transposase'
                transposase_sizes.append(size)  # Append the size to transposase_sizes list
        elif feature_type == "tRNA":  # If the feature is a tRNA
            trna_sizes.append(size)  # Append the size to trna_sizes list
        
        if intergenic_start is not None:  # If intergenic_start is not None
            intergenic_size = start - intergenic_start - 1  # Calculate the size of the intergenic region
            if intergenic_size > 0:  # If the size is positive
                intergenic_sizes.append(intergenic_size)  # Append the size to intergenic_sizes list
            intergenic_start = end  # Update intergenic_start to the end of the current feature
        else:
            intergenic_start = end  # Set intergenic_start to the end of the first feature

    return cds_sizes, trna_sizes, transposase_sizes, intergenic_sizes  # Return the lists of feature sizes

def analyze_mag(mag_id, fasta_folder, gff_folder):
    fasta_file = os.path.join(fasta_folder, f"{mag_id}.fna")  # Construct the path to the FASTA file
    gff_file = os.path.join(gff_folder, f"{mag_id}.gff")  # Construct the path to the GFF file

    total_assembly_size = calculate_total_assembly_size(fasta_file)  # Calculate the total assembly size using the FASTA file
    cds_sizes, trna_sizes, transposase_sizes, intergenic_sizes = process_gff(gff_file)  # Parse the GFF file to get sizes of features

    analysis_result = {  # Create a dictionary with the analysis results
        "MAG_ID": mag_id,  # MAG ID
        "Total_Assembly_Size_bp": total_assembly_size,  # Total assembly size in base pairs
        "CDS_Count": len(cds_sizes),  # Count of CDS features
        "CDS_Total_Size_bp": sum(cds_sizes),  # Total size of CDS features
        "CDS_Median_Size_bp": median(cds_sizes) if cds_sizes else 0,  # Median size of CDS features
        "tRNA_Count": len(trna_sizes),  # Count of tRNA features
        "tRNA_Total_Size_bp": sum(trna_sizes),  # Total size of tRNA features
        "tRNA_Median_Size_bp": median(trna_sizes) if trna_sizes else 0,  # Median size of tRNA features
        "Transposase_Count": len(transposase_sizes),  # Count of transposase features
        "Transposase_Total_Size_bp": sum(transposase_sizes),  # Total size of transposase features
        "Transposase_Median_Size_bp": median(transposase_sizes) if transposase_sizes else 0,  # Median size of transposase features
        "Intergenic_Count": len(intergenic_sizes),  # Count of intergenic regions
        "Intergenic_Total_Size_bp": sum(intergenic_sizes),  # Total size of intergenic regions
        "Intergenic_Median_Size_bp": median(intergenic_sizes) if intergenic_sizes else 0  # Median size of intergenic regions
    }

    return analysis_result  # Return the analysis result

def analyze_all_mags(fasta_folder, gff_folder):
    results = []  # List to store results for all MAGs
    for file_name in os.listdir(fasta_folder):  # Iterate over all files in the FASTA directory
        if file_name.endswith('.fa') or file_name.endswith('.fna'):  # Process files with extensions '.fa' or '.fna'
            mag_id = os.path.splitext(file_name)[0]  # Extract the MAG ID from the file name (without extension)
            try:
                result = analyze_mag(mag_id, fasta_folder, gff_folder)  # Analyze the MAG
                results.append(result)  # Append the result to the results list
            except Exception as e:  # Handle any exceptions that occur during the analysis
                print(f"Error processing {mag_id}: {e}")  # Print an error message

    return results  # Return the list of results

# Define the folders containing the FASTA and GFF files
fasta_folder = 'fasta'
gff_folder = 'gff'

# Analyze all MAGs in the specified folders
all_mag_results = analyze_all_mags(fasta_folder, gff_folder)

# Save results to a TSV file
df = pd.DataFrame(all_mag_results)  # Convert the results to a pandas DataFrame
output_file = 'all_mag_analysis_results.tsv'  # Define the output file name
df.to_csv(output_file, sep='\t', index=False)  # Save the DataFrame to a TSV file

print(f"Results saved to {output_file}")  # Print a message indicating that the results have been saved

