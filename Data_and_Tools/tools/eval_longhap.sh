laneid=$1
if [ -z $1 ]; then
    echo "please provide a library id"; exit 1
fi

vcf=${laneid}_sentieon.filtered.snp.vcf

perl /home/ysun/validate.pl $vcf ./LibOutDir | tee longhap_results.txt
perl /home/ysun/cal_ori_check.pl $vcf ./LibOutDir | tee -a longhap_results.txt
perl /home/ysun/evaluate.pl /home/ysun/GIAB_VCF/ ./LibOutDir | tee -a longhap_results.txt
