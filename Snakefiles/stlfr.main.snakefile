rule map_reads:
    input:
        ref = REF,
        fqs = expand("{sample}", sample=config['samples']['fastq'])
    output:
        bam = "Align/{}.sort.bam".format(config['samples']['id']), # may want this to be a temp file
        sam = "Align/{}.aln_mem.sam".format(config['samples']['id'])
    threads:
        config['threads']['bwa']
    params:
        sen_install = config['params']['sentieon_install'],
        sen_license = config['params']['sentieon_license'],
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['id'],
                                                          config['params']['platform'])
    benchmark:
        "Benchmarks/main.map_reads.txt"
    shell:
        "{params.sen_install}/bin/bwa mem -M -R '{params.readgroup}' -C "
            "-t {threads} {input.ref} {input.fqs} 2>Align/aln.err | tee {output.sam} | "
        "{params.sen_install}/bin/sentieon util sort -o {output.bam} -t {threads} --sam2bam -i -"


rule locus_collector:
    input:
        "Align/{id}.sort.bam"
    output:
        "Align/{id}_score.txt"
    threads:
        config['threads']['bwa']
    params:
        sen_install = config['params']['sentieon_install']
    benchmark:
        "Benchmarks/main.locus_collectord.{id}.txt"
    shell:
        "{params.sen_install}/bin/sentieon driver  -t {threads} "
            "-i {input} "
            "--algo LocusCollector "
            "--fun score_info {output}"


rule mark_dups:
    input:
        bam = "Align/{id}.sort.bam",
        score = "Align/{id}_score.txt"
    output:
        bam = "Align/{id}.sort.rmdup.bam",
        metrics = "Align/{id}_dedup_metrics.txt"
    threads:
        config['threads']['bwa']
    params:
        sen_install = config['params']['sentieon_install']
    benchmark:
        "Benchmarks/main.mark_dups.{id}.txt"
    shell:
        "{params.sen_install}/bin/sentieon driver  -t {threads} "
            "-i {input.bam} "
            "--algo Dedup "
            "--score_info {input.score} "
            "--metrics {output.metrics} {output.bam}"


rule mark_dups_txt:
    input:
        "Align/{id}_dedup_metrics.txt"
    output:
        "Align/{id}_dedup_metrics2.txt"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/main.mark_dups_txt.{id}.txt"
    shell:
        "perl {params.toolsdir}/tools/mark_dups_txt.pl {input} {output}"


rule run_longhap:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{}_sort.rmdup.addBX.bam".format(config['samples']['id']),
        vcf = "Make_Vcf/step1_haplotyper/{}_sentieon_pass_vars.vcf".format(config['samples']['id'])
    output:
        "Make_Vcf/step4_longhap/LibOutDir/chr1.vcf"
    threads:
        config['threads']['gnu_parallel']
    params:
        toolsdir = config['params']['toolsdir'],
        samp = config['samples']['id']
    benchmark:
        "Benchmarks/main.run_longhap.txt"
    shell:
        "mkdir -p Make_Vcf/step4_longhap; "
        "cd Make_Vcf/step4_longhap; "
        "sh {params.toolsdir}/tools/runlonghap.sh {params.samp} {threads}"


rule eval_longhap:
    input:
        "Make_Vcf/step4_longhap/LibOutDir/chr1.vcf"
    output:
        "Make_Vcf/step4_longhap/longhap_results.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        samp = config['samples']['id']
    benchmark:
        "Benchmarks/main.eval_longhap.txt"
    shell:
        "cd Make_Vcf/step4_longhap; "
        "sh {params.toolsdir}/tools/eval_longhap.sh {params.samp}"


#rule run_hapcut:
#    input:
#        sam = "Calc_Frag_Length/step1_removedup_rm000/{id}.sort.removedup_rm000.sam",
#        vcf = "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf"
#    output:
#        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/hapblock_{id}_chr1",
#        "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup.addBX.bam"
#    threads:
#        config['threads']['gnu_parallel']
#    params:
#        toolsdir = config['params']['toolsdir'],
#        samp = config['samples']['id']
#    shell:
#        "mkdir -p Make_Vcf/step3_hapcut; "
#        "cd Make_Vcf/step3_hapcut; "
#        "sh {params.toolsdir}/tools/run_hapcut2.sh {params.samp} {threads}"


rule run_hapcut_bam:
    input:
        "Calc_Frag_Length/step1_removedup_rm000/{id}.sort.removedup_rm000.sam",
    output:
        "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup.addBX.bam"
    threads:
        config['threads']['gnu_parallel']
    params:
        toolsdir = config['params']['toolsdir'],
        samp = config['samples']['id']
    benchmark:
        "Benchmarks/main.run_hapcut_bam.{id}.txt"
    shell:
        "mkdir -p Make_Vcf/step3_hapcut; "
        "cd Make_Vcf/step3_hapcut; "
        "sh {params.toolsdir}/tools/run_hapcut2_bam.sh {params.samp} {threads}"


rule run_hapcut_vcf:
    input:
        vcf = "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf",
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup.addBX.bam"
    output: 
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/hapblock_{id}_chr1"
    threads:
        config['threads']['gnu_parallel']
    params:
        toolsdir = config['params']['toolsdir'],
        samp = config['samples']['id']
    benchmark:
        "Benchmarks/main.run_hapcut_vcf.{id}.txt"
    shell:
        "mkdir -p Make_Vcf/step3_hapcut; "
        "cd Make_Vcf/step3_hapcut; "
        "sh {params.toolsdir}/tools/run_hapcut2_vcf.sh {params.samp} {threads}"

rule eval_hapcut:        
    input:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/hapblock_{}_chr1".format(config['samples']['id'])
    output:
        "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_comparison_with_giab.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        samp = config['samples']['id']
    benchmark:
        "Benchmarks/main.eval_hapcut.txt"
    shell:
        "mkdir -p Make_Vcf/step3_hapcut/step4_compare_with_refphasing; "
        "cd Make_Vcf/step3_hapcut; "
        "sh {params.toolsdir}/tools/compare_haps.sh {params.samp}"

        
