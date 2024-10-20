#!/bin/bash

# Path to your text file with SRR IDs
SRR_LIST="/home/vibhor/haseena_ip/mesoderm/SRR_Acc_List.txt"
# Read each SRR ID from the file and download the corresponding FASTQ.gz file
while IFS= read -r SRR_ID; do
  echo "Downloading and compressing $SRR_ID"
  fastq-dump --split-files --gzip "$SRR_ID"
done < "$SRR_LIST"

echo "Download complete."
