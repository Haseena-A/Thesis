import pandas as pd

# Define file paths
gene_list_file = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/genes_to_select_pancreas_cancer.txt'
bed_file = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/extracted_column_from_gtf_data.bed'
filtered_bed_file = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/filtered_pan_genes_with_chr.bed'

# Load the gene list from the text file
with open(gene_list_file, 'r') as file:
    gene_list = [line.strip().strip('"') for line in file]

# Load the BED file into a DataFrame
bed_df = pd.read_csv(bed_file, sep='\t', header=None)

# Filter the DataFrame based on the gene list
filtered_df = bed_df[bed_df[3].astype(str).isin(gene_list)]  # Ensure gene_list is string type for comparison

# Save the filtered DataFrame to a new BED file
filtered_df.to_csv(filtered_bed_file, sep='\t', header=False, index=False, float_format='%.0f')