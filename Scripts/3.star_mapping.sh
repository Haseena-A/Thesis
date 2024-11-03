for fastq_file in /data1/haseena_IP/tsv_files/sandisk/sandisk_output3/*.fastq.gz; do
    sample_name=$(basename $fastq_file .fastq.gz)

    STAR --runThreadN 8 \
         --genomeDir /data2/haseena/haseena/IP/hg38_index \
         --readFilesIn $fastq_file \
         --readFilesCommand zcat \
         --outFileNamePrefix /data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/${sample_name}_ \
         --outSAMtype BAM SortedByCoordinate
done
