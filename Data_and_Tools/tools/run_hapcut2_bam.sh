set -e
export laneid=$1
if [ -z $2  ]; then
    threads=10
else
    threads=$2
fi

align_sam=../../../Calc_Frag_Length/step1_removedup_rm000/${laneid}.sort.removedup_rm000.sam



###################################################################################################
###################################################################################################
linkdist=100000
chroms=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX
linkfragment_py=/research/rv-02/home/qmao/singletubesam/hapcut2_test/HapCUT2-master_20170728/HapCUT2-master/utilities/LinkFragments.py
####################################### vcf file dir ############################################
vcf_dir=../../../Make_Vcf/step1_haplotyper
vcf_prefix=${laneid}_sentieon_pass_vars
###################################################################################################


###step1: modify bam file to add barcode info as BX field
echo "begin step1"
mkdir -p step1_modify_bam
cd step1_modify_bam


awk -F $'\t' '{if($1~/#/){barcode=split($1,a,"#")
    print $0,"BX:Z:"a[2]}
    else{print}
    }' OFS="\t" $align_sam | samtools view -bh - > ${laneid}_sort.rmdup.addBX.bam


samtools index ${laneid}_sort.rmdup.addBX.bam
echo -e ${chroms//,/'\n'} | parallel -j $threads samtools view -bh ${laneid}_sort.rmdup.addBX.bam {} \> ${laneid}_sort.rmdup.addBX_{}.bam \&\& samtools index ${laneid}_sort.rmdup.addBX_{}.bam


cd ../
echo "step1 done"
date
