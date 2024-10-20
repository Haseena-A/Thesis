import pandas as pd

# Define the file path
input_file_path = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/filtered_pan_genes_with_chr.bed'
output_file_path = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/processed_pan_genes_with_chr.bed'

# Read the data from the file
df = pd.read_csv(input_file_path, sep="\t", header=None, names=['chr', 'start', 'end', 'name'])

# Group by the 4th column and find min and max of the 2nd and 3rd columns
result = df.groupby('name').agg({'chr': 'first', 'start': 'min', 'end': 'max'}).reset_index()

# Reorder columns to match BED format
result = result[['chr', 'start', 'end', 'name']]

# Save the result to a new file in BED format
result.to_csv(output_file_path, sep='\t', header=False, index=False)

print(f"Processed file saved to {output_file_path}")