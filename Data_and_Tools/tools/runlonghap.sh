#!/bin/bash
set -eo pipefail
script_name=`basename -s .sh "$0"`

laneid=$1

if [ -z $2  ]; then
    threads=10
else
    threads=$2
fi 


chroms=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX
bam=../../Align/${laneid}.sort.rmdup.bam
prevcf=../step1_haplotyper/${laneid}_sentieon_pass_vars.vcf
vcf=${laneid}_sentieon.filtered.snp.vcf

/opt/cgi-tools/bin/bcftools view -O v --type snps $prevcf > $vcf
logfile="${script_name}_$(date +"%Y%m%d%H%M").log"

echo "Generating LongHap Scripts"
perl /home/eanderson/LongHap_v1.3/Main.pl $bam $vcf 10000000 300000 1 1 LibShDir LibTmpDir LibOutDir

echo "Executing LongHap Scripts"
cd LibShDir
echo -e ${chroms//,/'\n'} | parallel -j $threads 'sh run.{}.sh'
cd ../

echo "LongHap Finished"

