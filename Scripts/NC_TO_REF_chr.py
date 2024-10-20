import re

# Function to convert RefSeq ID to "chr" format
def convert_chromosome_name(chrom_name):
    refseq_to_chr = {
        "NC_000001.11": "chr1",
        "NC_000002.12": "chr2",
        "NC_000003.12": "chr3",
        "NC_000004.12": "chr4",
        "NC_000005.10": "chr5",
        "NC_000006.12": "chr6",
        "NC_000007.14": "chr7",
        "NC_000008.11": "chr8",
        "NC_000009.12": "chr9",
        "NC_000010.11": "chr10",
        "NC_000011.10": "chr11",
        "NC_000012.12": "chr12",
        "NC_000013.11": "chr13",
        "NC_000014.9": "chr14",
        "NC_000015.10": "chr15",
        "NC_000016.10": "chr16",
        "NC_000017.11": "chr17",
        "NC_000018.10": "chr18",
        "NC_000019.10": "chr19",
        "NC_000020.11": "chr20",
        "NC_000021.9": "chr21",
        "NC_000022.11": "chr22",
        "NC_000023.11": "chrX",
        "NC_000024.10": "chrY",
    }
    return refseq_to_chr.get(chrom_name, chrom_name)

# Function to process the VCF file and convert chromosome names, writing to a new file
def process_and_save_vcf(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                columns = line.split('\t')
                columns[0] = convert_chromosome_name(columns[0])
                outfile.write('\t'.join(columns))

# Paths to the input and output VCF files
input_vcf = '/data1/haseena_IP/PRJNA431985_Breast_cancer/xenograft_ctc/ref1/All_xeno_merged_bam.vcf'
output_vcf ='/data1/haseena_IP/PRJNA431985_Breast_cancer/xenograft_ctc/ref1/All_xeno_merged_converted_bam.vcf'

# Process the input VCF file and save the converted version to a new file
process_and_save_vcf(input_vcf, output_vcf)

