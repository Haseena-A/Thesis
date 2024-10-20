input_file = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/processed_pan_genes_with_chr.bed'
output_file = '/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/processed_genes_without_last_column.bed'

try:
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Split the line by tab and take the first three columns
            columns = line.strip().split()
            if len(columns) >= 3:
                new_line = '\t'.join(columns[:3]) + '\n'
                # Write the new line to the output file
                outfile.write(new_line)
            else:
                print(f"Skipping line due to insufficient columns: {line}")
    print("Processed BED file saved successfully.")
except Exception as e:
    print(f"An error occurred: {e}")