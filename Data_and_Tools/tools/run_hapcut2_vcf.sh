set -eo pipefail
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

mkdir -p step2_split_vcf
cd step2_split_vcf


/opt/cgi-tools/bin/bcftools view -O z ${vcf_dir}/${vcf_prefix}.vcf > ${vcf_prefix}.vcf.gz
tabix -p vcf -f ${vcf_prefix}.vcf.gz
echo -e ${chroms//,/'\n'} | parallel -j $threads /opt/cgi-tools/bin/bcftools view -O v -r {} ${vcf_prefix}.vcf.gz \> ${vcf_prefix}_{}.vcf


cd ../
date

vcf_dir=../step2_split_vcf


##use hapcut2 as of 20170731 version
echo "begin step2"
mkdir -p step3_run_hapcut2_10xpipeline
cd step3_run_hapcut2_10xpipeline


mkdir -p s1_unlinked_frag
mkdir -p s2_link_frag_files
mkdir -p s3_hapcut_output

for i in {chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX}
do
echo $i

bamfile=../step1_modify_bam/${laneid}_sort.rmdup.addBX_${i}.bam

echo "date && extractHAIRS --10X 1 --bam "$bamfile" --VCF "${vcf_dir}"/"${vcf_prefix}"_"${i}".vcf --out s1_unlinked_frag/unlinked_fragment_"${laneid}"_"${i}" && python3 "$linkfragment_py" --bam "$bamfile" --vcf "${vcf_dir}"/"${vcf_prefix}"_"${i}".vcf --fragments s1_unlinked_frag/unlinked_fragment_"${laneid}"_"${i}" --out s2_link_frag_files/linked_frag_"${laneid}"_"${i}" -d "$linkdist" && HAPCUT2 --nf 1 --fragments s2_link_frag_files/linked_frag_"${laneid}"_"${i}" --vcf "${vcf_dir}"/"${vcf_prefix}"_"${i}".vcf --output s3_hapcut_output/hapblock_"${laneid}"_"${i}" && date" > run_hapcut_for_${i}.sh

done

echo -e ${chroms//,/'\n'} | parallel -j $threads nohup sh run_hapcut_for_{}.sh \> hapcut_${laneid}_{}.log

wait
cd ../
echo "step2 done"
date
