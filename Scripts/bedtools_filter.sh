# Define the directory paths
input_dir="/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM"
output_dir="/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/filtered_bam"
bed_file="/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/processed_genes_without_last_column.bed"

# Loop through each BAM file in the input directory
for bam_file in "$input_dir"/*.bam; do
    # Extract the SRR ID from the BAM file name
    srr_id=$(basename "$bam_file" | cut -d'_' -f1)
    
    # Define the output BAM file name
    output_file="$output_dir/${srr_id}_filtered.bam"
    
    # Run bedtools intersect
    bedtools intersect -abam "$bam_file" -b "$bed_file" > "$output_file"
    
    echo "Processed $bam_file and saved filtered output to $output_file"
done

echo "All BAM files processed successfully."