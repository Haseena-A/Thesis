#!/bin/bash

GTF=/home/vibhor/haseena_ip/hg19_files/hg19.ncbiRefSeq.gtf
FILES=/home/vibhor/haseena_ip/stem_bam/BAM1/*.bam
OUTPUT=/home/vibhor/haseena_ip/stem_bam/loom

for f in $FILES
do
   echo $f
   velocyto run -c -U -o $OUTPUT -m /home/vibhor/haseena_ip/hg19_files/hg19_rmsk.gtf $f $GTF
done

echo "done!"
