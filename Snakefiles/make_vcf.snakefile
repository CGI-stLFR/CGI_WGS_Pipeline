
rule run_dnascope:
    input:
        bam = "Align/{id}.sort.rmdup.bam",
        ref = REF
    output:
        "Make_Vcf/step1_haplotyper/{id}_sentieon.tmp.vcf"
    threads:
        config['threads']['haplotyper']
    params:
        sen_install = config['params']['sentieon_install'],
        sen_model = config['params']['sentieon_model'],
        dbsnp = config['params']['dbsnp_path']
    shell:
        "{params.sen_install}/bin/sentieon driver -r {input.ref}  -t {threads} "
            "-i {input.bam} "
            "--algo DNAscope "
            "-d {params.dbsnp} "
            "--model {params.sen_model} "            
            "{output}"


rule run_dnamodel_apply:
    input:
        vcf = "Make_Vcf/step1_haplotyper/{id}_sentieon.tmp.vcf",
        ref = REF
    output:
        "Make_Vcf/step1_haplotyper/{id}_sentieon.vcf"
    threads:
        config['threads']['haplotyper']
    params:
        sen_install = config['params']['sentieon_install'],
        sen_model = config['params']['sentieon_model']
    shell:
        "{params.sen_install}/bin/sentieon driver -t {threads} -r {input.ref} "
        "--algo DNAModelApply "
        "--model {params.sen_model} "
        "-v {input.vcf} {output}"


rule keep_pass_vars:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon.vcf"
    output:
        "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf"
    shell:
        """
        awk '($1~/^#/ || $7=="PASS" || $7=="."){{print}}' {input} > {output}
        """

rule select_snps:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon.vcf"
    output:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz"
    shell:
        "/opt/cgi-tools/bin/bcftools view -O z --type snps {input} > {output}"


rule select_indels:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon.vcf"
    output:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz"
    shell:
        "/opt/cgi-tools/bin/bcftools view -O z --type indels {input} > {output}"


rule index_snps:
    input:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz"
    output:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz.tbi"
    shell:
        "/opt/cgi-tools/bin/tabix -p vcf -f {input}"


rule index_indels:
    input:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz"
    output:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz.tbi"
    shell:
        "/opt/cgi-tools/bin/tabix -p vcf -f {input}"

rule eval_snps:
    input:
        snp = "Make_Vcf/step2_benchmarking/{}.snp.vcf.gz".format(config['samples']['id']),
        index = "Make_Vcf/step2_benchmarking/{}.snp.vcf.gz.tbi".format(config['samples']['id'])
    output:
        "Make_Vcf/step2_benchmarking/snp_compare/summary.txt"
    params:
        benchmark_snp = config['benchmark']['benchmark_snp'],
        ref_sdf = config['benchmark']['ref_sdf'],
        bedfile = config['benchmark']['bedfile'],
        direc = "Make_Vcf/step2_benchmarking/snp_compare"
    shell:
        "rm -r {params.direc};"
        "/research/rv-02/home/qmao/Scripts/rtg-tools-3.8.4/rtg vcfeval "
            "-b {params.benchmark_snp} "
            "-c {input.snp} "
            "-e {params.bedfile} "
            "-t {params.ref_sdf} "
            "-o {params.direc}"


rule eval_indels:
    input:
        indel = "Make_Vcf/step2_benchmarking/{}.indel.vcf.gz".format(config['samples']['id']),
        index = "Make_Vcf/step2_benchmarking/{}.indel.vcf.gz.tbi".format(config['samples']['id'])
    output:
        "Make_Vcf/step2_benchmarking/indel_compare/summary.txt"
    params:
        benchmark_indel = config['benchmark']['benchmark_indel'],
        ref_sdf = config['benchmark']['ref_sdf'],
        bedfile = config['benchmark']['bedfile'],
        direc = "Make_Vcf/step2_benchmarking/indel_compare"
    shell:
        "rm -r {params.direc};"
        "/research/rv-02/home/qmao/Scripts/rtg-tools-3.8.4/rtg vcfeval "
            "-b {params.benchmark_indel} "
            "-c {input.indel} "
            "-e {params.bedfile} "
            "-t {params.ref_sdf} "
            "-o {params.direc} "

