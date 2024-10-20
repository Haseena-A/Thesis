import os
import gzip

def load_barcodes(barcode_file):
    """Load the valid spot-specific barcodes from a file."""
    with open(barcode_file, 'r') as f:
        barcodes = {line.strip() for line in f}
    return barcodes

def demultiplex_fastq(r1_path, r2_path, barcodes, output_dir):
    """Demultiplex the gzipped FASTQ files based on the spot-specific barcodes."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with gzip.open(r1_path, 'rt') as r1, gzip.open(r2_path, 'rt') as r2:
        while True:
            # Read the next 4 lines for each FASTQ read in R1 and R2
            r1_lines = [r1.readline().strip() for _ in range(4)]
            r2_lines = [r2.readline().strip() for _ in range(4)]
            
            # If end of file is reached
            if not r1_lines[0] or not r2_lines[0]:
                break

            # Extract barcode and UMI from R1
            barcode = r1_lines[1][:18]  # First 18 nt for the spot-specific barcode
            umi = r1_lines[1][18:25]  # Next 7 nt for the UMI

            # Check if the barcode is valid
            if barcode in barcodes:
                # Tag R2 read with barcode and UMI
                r2_tagged_header = f"{r2_lines[0]}|{barcode}|{umi}"
                
                # Prepare the R2 read with the tagged header
                r2_tagged = [r2_tagged_header] + r2_lines[1:]

                # Create the output file if it doesn't exist
                output_file = os.path.join(output_dir, f"{barcode}_R2.fastq.gz")
                with gzip.open(output_file, 'at') as out_f:
                    for line in r2_tagged:
                        out_f.write(line + '\n')

def main():
    barcode_file = "/data1/haseena_IP/spatial_data/star_output/GSE111672_barcodes_only.txt"  # File containing valid barcodes, each 18 nt long
    r1_fastq = "/data1/haseena_IP/spatial_data/trim/trim_galore_output/SRR10197468_1_val_1_val_1.fq.gz"  # Input R1 FASTQ file (gzipped)
    r2_fastq = "/data1/haseena_IP/spatial_data/trim/trim_galore_output/SRR10197468_2_val_2_val_2.fq.gz"  # Input R2 FASTQ file (gzipped)
    output_directory = "/data1/haseena_IP/tsv_files/sandisk/sandisk_output3"  # Directory to store demultiplexed FASTQs

    # Load valid barcodes
    barcodes = load_barcodes(barcode_file)
    
    # Run demultiplexing
    demultiplex_fastq(r1_fastq, r2_fastq, barcodes, output_directory)

if __name__ == "__main__":
    main()