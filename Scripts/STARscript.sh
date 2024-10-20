
for fastq_file_1 in /home/vibhor/haseena_ip/mesoderm/*_1.fastq.gz; do
    sample_name=$(basename $fastq_file_1 _1.fastq.gz)
    fastq_file_2=/home/vibhor/haseena_ip/mesoderm/${sample_name}_2.fastq.gz

    STAR --runThreadN 8 \
         --genomeDir /home/vibhor/haseena_ip/Genome_index_hg19 \
         --readFilesIn $fastq_file_1 $fastq_file_2 \
         --readFilesCommand zcat \
         --outFileNamePrefix /home/vibhor/haseena_ip/mesoderm/BAM/${sample_name}_ \
         --outSAMtype BAM SortedByCoordinate
done