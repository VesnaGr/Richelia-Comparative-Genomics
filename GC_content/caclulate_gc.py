# This script processes genome files and their corresponding GFF annotation files to extract coding sequences (CDS) and intergenic sequences (IGS).
# It then computes the GC content for these sequences and saves the results in separate TSV files for each genome.

from Bio import SeqIO # Used for reading and writing sequence files
from Bio.SeqRecord import SeqRecord # Used for creating sequence records
from Bio.SeqUtils import GC # Used for computing the GC content of sequences
from BCBio import GFF # Used for parsing GFF files
import os # Used for file and directory operations

# This function extracts CDS and IGS regions from a genome sequence based on annotations in a GFF file and writes them to separate FASTA files
def extract_regions_from_genome(genome_record, gff_path, output_cds, output_intergenic, genome_records_dict):
    cds_records = [] # Initializes an empty list for CDS records
    intergenic_records = [] # Initializes an empty list for intergenic records
    prev_feature_end = 0 # Sets the end position of the previous feature to 0

    print(f"Processing GFF file: {gff_path}")
    with open(gff_path, "r") as gff_file: # Opens the file specified by 'gff_path' in read mode ("r")
        for record in GFF.parse(gff_file, base_dict=genome_records_dict): # Parses the GFF file
            seq_id = record.id
            genome_record = genome_records_dict[seq_id]

            for feature in sorted(record.features, key=lambda x: x.location.start): # Sorts features by their start location
                if feature.type == "gene":
                    # Add the region from the end of the last gene to the start of this one as intergenic
                    if prev_feature_end < feature.location.start:
                        intergenic_records.append(SeqRecord(genome_record.seq[prev_feature_end:feature.location.start], id=f"intergenic_{prev_feature_end}_{feature.location.start}", description=""))

                    # Add the gene to the CDS records
                    cds_records.append(SeqRecord(genome_record.seq[feature.location.start:feature.location.end], id=feature.id, description=""))

                    # Update the previous feature end
                    prev_feature_end = feature.location.end

            # Capture the region from the end of the last gene to the end of the genome record
            if prev_feature_end < len(genome_record.seq):
                intergenic_records.append(SeqRecord(genome_record.seq[prev_feature_end:], id=f"intergenic_{prev_feature_end}_{len(genome_record.seq)}", description=""))

    # Write records to files
    print(f"Writing CDS to {output_cds}")
    SeqIO.write(cds_records, output_cds, "fasta") # Writes the CDS records to a FASTA file
    print(f"Writing intergenic regions to {output_intergenic}")
    SeqIO.write(intergenic_records, output_intergenic, "fasta") # Writes the intergenic records to a FASTA file

# This function computes the GC content for short and long intergenic sequences and CDS sequences, and saves the results in a TSV file
def compute_gc_content(intergenic_path, cds_path, output_gc):
    short_intergenic_gc = [] # Initializes an empty list for short intergenic GC content
    long_intergenic_gc = [] # Initializes an empty list for long intergenic GC content
    cds_gc = [] # Initializes an empty list for CDS GC content

    # Compute GC content for intergenic regions
    print(f"Computing GC content for intergenic regions from {intergenic_path}")
    for record in SeqIO.parse(intergenic_path, "fasta"):
        gc_content = GC(record.seq) # Computes the GC content of the sequence
        if len(record.seq) < 300:
            short_intergenic_gc.append((record.id, gc_content)) # Adds the GC content to the short intergenic list if the sequence length is less than 300 base pairs
        else:
            long_intergenic_gc.append((record.id, gc_content)) # Adds the GC content to the long intergenic list if the sequence length is 300 base pairs or more

    # Compute GC content for CDS regions
    print(f"Computing GC content for CDS regions from {cds_path}")
    for record in SeqIO.parse(cds_path, "fasta"):
        gc_content = GC(record.seq) # Computes the GC content of the sequence
        cds_gc.append((record.id, gc_content)) # Adds the GC content to the CDS list

    # Saving results to a TSV file
    print(f"Saving GC content results to {output_gc}")
    with open(output_gc, "w") as output_file:
        # Write header
        output_file.write("Category\tSeq_ID\tGC_Content\n")

        # Short intergenic sequences
        for seq_id, gc_value in short_intergenic_gc:
            output_file.write(f"Short Intergenic\t{seq_id}\t{gc_value:.2f}\n")

        # Long intergenic sequences
        for seq_id, gc_value in long_intergenic_gc:
            output_file.write(f"Long Intergenic\t{seq_id}\t{gc_value:.2f}\n")

        # CDS sequences
        for seq_id, gc_value in cds_gc:
            output_file.write(f"CDS\t{seq_id}\t{gc_value:.2f}\n")

# Define the directories
genome_folder = "fasta" # Directory containing genome FASTA files
gff_folder = "gff" # Directory containing GFF files
output_folder = "output" # Directory to save the output files

# Ensure the output directory exists
os.makedirs(output_folder, exist_ok=True) # Creates the output directory if it does not exist

# Process each genome file in the folder
for genome_filename in os.listdir(genome_folder):
    if genome_filename.endswith('.fa') or genome_filename.endswith('.fna'):
        genome_path = os.path.join(genome_folder, genome_filename)
        genome_id = os.path.splitext(genome_filename)[0]
        
        gff_path = os.path.join(gff_folder, f"{genome_id}.gff")
        output_cds_path = os.path.join(output_folder, f"{genome_id}_cds.fasta")
        output_intergenic_path = os.path.join(output_folder, f"{genome_id}_intergenic.fasta")
        output_gc_path = os.path.join(output_folder, f"{genome_id}_gc_content.tsv")
        
        # Read the genome sequence
        print(f"Reading genome file: {genome_path}")
        genome_records_dict = {rec.id: rec for rec in SeqIO.parse(genome_path, "fasta")} # Reads the genome sequence and stores it in a dictionary
        
        # Process each genome record and its annotations
        for genome_record_id, genome_record in genome_records_dict.items():
            extract_regions_from_genome(genome_record, gff_path, output_cds_path, output_intergenic_path, genome_records_dict)
        
        # Compute and save GC content
        compute_gc_content(output_intergenic_path, output_cds_path, output_gc_path)

print("Processing complete.")

