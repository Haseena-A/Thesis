fastq_file=/data1/haseena_IP/prostate_cancer_fastq/sample.fastq.gz
sample_name=$(basename "$fastq_file" .fastq.gz)

STAR --runThreadN 8 \
     --genomeDir /data2/haseena/haseena/IP/hg38_index \
     --readFilesIn "$fastq_file" \
     --readFilesCommand zcat \
     --outFileNamePrefix /data1/haseena_IP/prostate_cancer_bam/"${sample_name}_" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --quantMode TranscriptomeSAM GeneCounts \
     --twopassMode Basic \
     --outFilterMismatchNmax 10 \
     --outFilterMismatchNoverLmax 0.1 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000
