wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar axvf annovar.latest.tar.gz
cd annovar

perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 refGene humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 refGeneWithVer humandb/
#perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg19 refGene humandb/
#perl annotate_variation.pl --buildver hg19 --downdb seq humandb/hg19_seq
#perl retrieve_seq_from_fasta.pl humandb/hg19_refGene.txt -seqdir humandb/hg19_seq -format refGene -outfile humandb/hg19_refGeneMrna.fa
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 1000g2015aug humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 esp6500siv2_all humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 esp6500siv2_aa humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 esp6500siv2_ea humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 exac03 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 nci60 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 avsnp147 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 cosmic70 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 clinvar_20220320 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 dbnsfp42a humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 gnomad_exome humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 gnomad_genome humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 dbscsnv11 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 ensGene humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 knownGene humandb/
#perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg19 ensGene humandb/
#perl retrieve_seq_from_fasta.pl humandb/hg19_ensGene.txt -seqdir humandb/hg19_seq -format ensGene -outfile humandb/hg19_ensGeneMrna.fa
#perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg19 knownGene humandb/
#perl retrieve_seq_from_fasta.pl humandb/hg19_knownGene.txt -seqdir humandb/hg19_seq -format knownGene -outfile humandb/hg19_knownGeneMrna.fa
perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg19 rmsk humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 regsnpintron humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 gene4denovo201907 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg19 icgc28 humandb/

perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 refGene humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 refGeneWithVer humandb/
#perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg38 refGene humandb/
#perl annotate_variation.pl --buildver hg38 --downdb seq humandb/hg38_seq
#perl retrieve_seq_from_fasta.pl humandb/hg38_refGene.txt -seqdir humandb/hg38_seq/chroms -format refGene -outfile humandb/hg38_refGeneMrna.fa
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 1000g2015aug humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 esp6500siv2_all humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 esp6500siv2_aa humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 esp6500siv2_ea humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 exac03 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 nci60 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 avsnp147 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 cosmic70 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 clinvar_20220320 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 dbnsfp42a humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 gnomad_exome humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 gnomad_genome humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 dbscsnv11 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 ensGene humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 knownGene humandb/
#perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg38 ensGene humandb/
#perl retrieve_seq_from_fasta.pl humandb/hg38_ensGene.txt -seqdir humandb/hg38_seq/chroms -format ensGene -outfile humandb/hg38_ensGeneMrna.fa
#perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg38 knownGene humandb/
#perl retrieve_seq_from_fasta.pl humandb/hg38_knownGene.txt -seqdir humandb/hg38_seq/chroms -format knownGene -outfile humandb/hg38_knownGeneMrna.fa
perl annotate_variation.pl --downdb --webfrom ucsc --buildver hg38 rmsk humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 gene4denovo201907 humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 regsnpintron humandb/
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 icgc28 humandb/

cd ..

perl annovar/table_annovar.pl annovar/example/ex2.vcf annovar/humandb/ --outfile ex2 --buildver hg38 \
    --protocol refGene,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,nci60,avsnp147,cosmic70,clinvar_20220320,dbnsfp42a,gnomad_exome,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene,regsnpintron,gene4denovo201907,icgc28 \
    --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,g,g,f,f,f \
    --vcfinput --otherinfo --thread $(nproc) --maxgenethread $(nproc)

git clone https://github.com/WGLab/InterVar.git -b v2.2.1

wget https://omim.org/static/omim/data/mim2gene.txt

mv mim2gene.txt InterVar/intervardb

python InterVar/Intervar.py \
    --input=ex2.avinput \
    --output=ex2 \
    --buildver=hg38 \
    --database_intervar=InterVar/intervardb \
    --table_annovar=annovar/table_annovar.pl \
    --convert2annovar=annovar/convert2annovar.pl \
    --annotate_variation=annovar/annotate_variation.pl \
    --database_locat=annovar/humandb \
    --skip_annovar