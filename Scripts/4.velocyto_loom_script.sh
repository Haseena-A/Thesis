#!/bin/bash

GTF=/home/haseena/hg19_files/hg19.ncbiRefSeq.gtf
FILES=/home/haseena/new_mesenchymal/0h_bam/demultiplexed_bam/*.bam
OUTPUT=/home/haseena/new_mesenchymal/0h_bam/demultiplexed_bam/loom

for f in $FILES
do
   echo $f
   velocyto run -c -U -o $OUTPUT -m /home/haseena/hg19_files/hg19_rmsk.gtf $f $GTF
done

echo "done!"
