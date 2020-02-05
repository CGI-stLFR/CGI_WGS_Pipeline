set -e

laneid=$1
if [ -z $1 ]; then
    echo "Please provide a library id"; exit 1
fi

###################################################################################################
###################################################################################################
linkdist=100000
chroms=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX
linkfragment_py=/research/rv-02/home/qmao/singletubesam/hapcut2_test/HapCUT2-master_20170728/HapCUT2-master/utilities/LinkFragments.py
contig_size=/research/rv-02/home/qmao/singletubesam/hapcut2_test/HapCUT2-master/HapCUT2-master/utilities/hg19.chrom.sizes
compare_block_py=/research/rv-02/home/qmao/singletubesam/hapcut2_test/HapCUT2-master_20170728/HapCUT2-master/utilities/calculate_haplotype_statistics.py

####################################### vcf file dir ############################################
vcf_dir=../step2_split_vcf
vcf_prefix=${laneid}_sentieon_pass_vars
align_sam=../../../Calc_Frag_Length/step1_removedup_rm000/${laneid}.sort.removedup_rm000.sam

################################phased vcf used as reference dir ###################################
ref_phased_vcf_dir=/research/rv-02/home/qmao/singletubesam/hapcut2_test/NA12878_vcf_giab
ref_phased_vcf_prefix=giab

##########################################
echo "begin step4"
mkdir -p step4_compare_with_refphasing
cd step4_compare_with_refphasing

output_file=hapcut_comparison_with_${ref_phased_vcf_prefix}.txt
h1_prefix=../step3_run_hapcut2_10xpipeline/s3_hapcut_output/hapblock_${laneid}
v1_prefix=${vcf_dir}/${vcf_prefix}
f1_prefix=../step3_run_hapcut2_10xpipeline/s2_link_frag_files/linked_frag_${laneid}

pv_prefix=${ref_phased_vcf_dir}/${ref_phased_vcf_prefix}

echo "compare ${laneid} with ${ref_phased_vcf_prefix}" > ${output_file}

for i in {chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX}
do

echo $i >> ${output_file}

python3 ${compare_block_py} -h1 ${h1_prefix}_${i} -v1 ${v1_prefix}_${i}.vcf -f1 ${f1_prefix}_${i} -pv ${pv_prefix}_${i}.vcf -c $contig_size >> ${output_file}

done

#######################
echo "combine all chrs" >> ${output_file}

python3 ${compare_block_py} -h1 ${h1_prefix}_chr1 ${h1_prefix}_chr2 ${h1_prefix}_chr3 ${h1_prefix}_chr4 ${h1_prefix}_chr5 ${h1_prefix}_chr6 ${h1_prefix}_chr7 ${h1_prefix}_chr8 ${h1_prefix}_chr9 ${h1_prefix}_chr10 ${h1_prefix}_chr11 ${h1_prefix}_chr12 ${h1_prefix}_chr13 ${h1_prefix}_chr14 ${h1_prefix}_chr15 ${h1_prefix}_chr16 ${h1_prefix}_chr17 ${h1_prefix}_chr18 ${h1_prefix}_chr19 ${h1_prefix}_chr20 ${h1_prefix}_chr21 ${h1_prefix}_chr22 ${h1_prefix}_chrX -v1 ${v1_prefix}_chr1.vcf ${v1_prefix}_chr2.vcf ${v1_prefix}_chr3.vcf ${v1_prefix}_chr4.vcf ${v1_prefix}_chr5.vcf ${v1_prefix}_chr6.vcf ${v1_prefix}_chr7.vcf ${v1_prefix}_chr8.vcf ${v1_prefix}_chr9.vcf ${v1_prefix}_chr10.vcf ${v1_prefix}_chr11.vcf ${v1_prefix}_chr12.vcf ${v1_prefix}_chr13.vcf ${v1_prefix}_chr14.vcf ${v1_prefix}_chr15.vcf ${v1_prefix}_chr16.vcf ${v1_prefix}_chr17.vcf ${v1_prefix}_chr18.vcf ${v1_prefix}_chr19.vcf ${v1_prefix}_chr20.vcf ${v1_prefix}_chr21.vcf ${v1_prefix}_chr22.vcf ${v1_prefix}_chrX.vcf -f1 ${f1_prefix}_chr1 ${f1_prefix}_chr2 ${f1_prefix}_chr3 ${f1_prefix}_chr4 ${f1_prefix}_chr5 ${f1_prefix}_chr6 ${f1_prefix}_chr7 ${f1_prefix}_chr8 ${f1_prefix}_chr9 ${f1_prefix}_chr10 ${f1_prefix}_chr11 ${f1_prefix}_chr12 ${f1_prefix}_chr13 ${f1_prefix}_chr14 ${f1_prefix}_chr15 ${f1_prefix}_chr16 ${f1_prefix}_chr17 ${f1_prefix}_chr18 ${f1_prefix}_chr19 ${f1_prefix}_chr20 ${f1_prefix}_chr21 ${f1_prefix}_chr22 ${f1_prefix}_chrX -pv ${pv_prefix}_chr1.vcf ${pv_prefix}_chr2.vcf ${pv_prefix}_chr3.vcf ${pv_prefix}_chr4.vcf ${pv_prefix}_chr5.vcf ${pv_prefix}_chr6.vcf ${pv_prefix}_chr7.vcf ${pv_prefix}_chr8.vcf ${pv_prefix}_chr9.vcf ${pv_prefix}_chr10.vcf ${pv_prefix}_chr11.vcf ${pv_prefix}_chr12.vcf ${pv_prefix}_chr13.vcf ${pv_prefix}_chr14.vcf ${pv_prefix}_chr15.vcf ${pv_prefix}_chr16.vcf ${pv_prefix}_chr17.vcf ${pv_prefix}_chr18.vcf ${pv_prefix}_chr19.vcf ${pv_prefix}_chr20.vcf ${pv_prefix}_chr21.vcf ${pv_prefix}_chr22.vcf ${pv_prefix}_chrX.vcf -c $contig_size >> ${output_file}


cd ../

echo "step4 done"
date
