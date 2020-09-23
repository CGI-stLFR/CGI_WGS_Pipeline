# generate bam for hapcut
rule make_hapcut_bam:
    input:
        "Align/{id}.sort.removedup_rm000.sam"
    output:
        "Align/{id}.sort.removedup_rm000.bam"
    shell:
        "samtools view -bh {input} > {output} && "
        "samtools index {output}"


# split the bam by chromosome
rule split_hapcut_bam:
    input:
        "Align/{id}.sort.removedup_rm000.bam"
    output:
        "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup_{chr}.bam"
    params:
        id = config['samples']['id']
    shell:
        "samtools view -bh {input} {wildcards.chr} > "
            "{output} && "
            "samtools index {output}"


# gzip the pass vars vcf
rule gzip_hapcut_vcf:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf"
    output:
        gz_vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz",
        index = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz.tbi"
    params:
        bcftools = config['params']['bcftools']
    shell:
        "{params.bcftools} view -O z {input} > "
            "{output.gz_vcf} && "
            "tabix -p vcf -f {output.gz_vcf}"


# split the gzipped vcf by chromosome
rule split_hapcut_vcf:
    input:
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz",
        index = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz.tbi"
    output:
        "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf"
    params:
        bcftools = config['params']['bcftools'],
        id = config['samples']['id']
    shell:
        "{params.bcftools} view -O v "
            "-r {wildcards.chr} "
            "{input.vcf} > "
            "{output}"


# run extract hairs to get fragments for hapcut
rule get_hapcut_fragments:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup_{chr}.bam",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf"
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    params:
        toolsdir = config['params']['toolsdir']
    shell:
        "{params.toolsdir}/tools/HapCUT2/build/extractHAIRS --10X 1 "
            "--bam {input.bam} "
            "--VCF {input.vcf} "
            "--out {output}"


# run LinkFragments.py
# This preps files for hapcut
rule link_hapcut_fragments:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup_{chr}.bam",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf",
        frag = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}"
    params:
        toolsdir = config['params']['toolsdir'],
        linkdist = config['params']['hapcut_link_dist']
    shell:
        "python3 {params.toolsdir}/tools/HapCUT2/utilities/LinkFragments.py --bam {input.bam} "
            "--VCF {input.vcf} "
            "--fragments {input.frag} "
            "--out {output} "
            "-d {params.linkdist}"


# Run hapcut for each chromosome
rule run_hapcut2:
    input:
        frag = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf"
    output:
        blocks = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}",
        vcf = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF"
    params:
        toolsdir = config['params']['toolsdir']
    shell:
        "{params.toolsdir}/tools/HapCUT2/build/HAPCUT2 --nf 1 "
            "--fragments {input.frag} "
            "--VCF {input.vcf} "
            "--output {output.blocks}"


# evaluate hapcut
rule evaluate_hapcut2:
    input:
        hapblock = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}", id=config['samples']['id'], chr=CHROMS),
        frag = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}", id=config['samples']['id'], chr=CHROMS),
        vcf = expand("Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf", id=config['samples']['id'], chr=CHROMS)
    output:
        "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        truth_chroms = config['benchmark']['truth_vcf_dir'],
        ref_index = REF + '.fai'
    run:
        import sys
        from pathlib import Path
        
        # if there's a truth_chroms parameter, we should be evaluating a benchmarked sample
        # the truth_chroms_{chr}.vcf should contain the phased truthset VCF for that chromosome
        if params.truth_chroms:
            truth_chroms = expand(params.truth_chroms + "/truth_chroms_{chr}.vcf", chr=CHROMS)
            shell("""python3 {toolsdir}/calculate_haplotype_statistics.py -h1 {input.hapblock} 
                     -v1 {input.vcf} -f1 {input.frag} 
                     -pv {truth_chroms} -c {params.ref_index} > {output}""")
        # Otherwise we import the hapblock_vcf_error_rate_multiple module from hapcut
        # and then use our vcf as the truth VCF
        # This means we won't get any switch statistics, but without knowing the truth 
        # we can't figure that out anyway
        # This does allow us to get the number of phased SNPs, N50 and NA50
        else:
            sys.path.append(str(Path(params.toolsdir) / 'tools'))
            from calculate_haplotype_statistics_no_truth import hapblock_vcf_error_rate_multiple

            err = hapblock_vcf_error_rate_multiple(input.hapblock, input.vcf, input.vcf, False)
            with open(output[0], "w") as outf:
                print(err, file=outf)


# Aggregate all of the hapcut VCFs as a full phased VCF
rule aggregate_hapcut_vcf:
    input:
        vcfs = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF", id=config['samples']['id'], chr=CHROMS)
    output:
        agg_vcf = "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id'])
    run:
        import sys, re, glob
        from pathlib import Path

        contigs = []
        # open up the first vcf to get the header and output it to the aggregte vcf
        with open(input.vcfs[0], "r") as vcf_header, open(output.agg_vcf, "w") as outfile:
            for line in vcf_header:
                if line.startswith("#"):
                    print(line.strip(), file=outfile)
                    # scrape the header for contigs present
                    if "contig=" in line:
                        # capture contig info
                        id_name = re.split("<|>", line.strip())[1]
                        # split into a dict of ID: contig name, length: <len>, assembly: ref
                        keys = dict(parts.split("=") for parts in id_name.split(","))
                        contigs.append(keys['ID'])
                else:
                    break

        
        file_dict = {}
        # iterate through input VCFs
        # to create a dictionary of contigs and file paths
        for infile in input.vcfs:
            name = Path(infile).name
            tail = name.split("_hapblock_")[-1]
            contig = tail.split(".phased.VCF")[0]
            file_dict[contig] = infile

        # iterate through contigs, this means the results will be in order
        for contig in contigs:
            # some contigs won't be present as files if they didn't have SNPs
            # this is more relevant for alignments to de novo assemblies
            try:
                infile = file_dict[contig]
            except KeyError:
                print(f"No file for contig {contig}, skipping...", file=sys.stderr)
                continue
            # if we find the file, open it and output results to the aggregate file
            with open(infile, "r") as vcf, open(output.agg_vcf, "a") as outfile:
                for line in vcf:
                    if line.startswith("#"):
                        continue
                    else:
                        print(line.strip(), file=outfile)


# get snps from the pass vars VCF for longhap
# hapcut won't try and phase indels but longhap will
# it also doesn't do great at phasing the indels
# so it's best not to have them at all
rule get_snps_for_longhap:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf"
    output:
        "Make_Vcf/step4_longhap/{id}_sentieon_pass_vars.vcf"
    params:
        bcftools = config['params']['bcftools']
    shell:
        "{params.bcftools} view -O v --type snps {input} > {output}"


# First we need to generate the scripts longhap uses to evaluate each chromosome
rule generate_longhap_scripts:
    input:
        bam = "Align/{}.sort.rmdup.bam".format(config['samples']['id']),
        vcf = "Make_Vcf/step4_longhap/{}_sentieon_pass_vars.vcf".format(config['samples']['id'])
    output:
        "Make_Vcf/step4_longhap/LibShDir/run.{chr}.sh"
    params:
        win_size = config['longhap']['win_size'],
        max_len = config['longhap']['max_lfr_len'],
        min_bc = config['longhap']['min_barcode'],
        min_link = config['longhap']['min_link'],
        longhap = config['params']['longhap']
    run:
        from pathlib import Path
        import os

        root = Path("Make_Vcf/step4_longhap")
        dirs = ["LibShDir", "LibTmpDir", "LibOutDir"]

        # create all the directories that LongHap expects
        for direc in dirs:
            direc_path = root/direc
            if not direc_path.exists():
                direc_path.mkdir()
        
        # define paths to the various longhap executables we need
        longhap_ex = Path(params.longhap)
        # the below script merges the various bins
        longhap_merge = longhap_ex.parent / "merge.pl"
        # this script creates a VCF from longhaps log files
        longhap_log2vcf = longhap_ex.parent / "log2vcf.pl"

        ref_index = REF + ".fai"
        ref_dict = {}
        
        # get the contig lengths from the fasta index
        # These are used to create the various bins for the shell scripts
        with open(ref_index, "r") as ref_ind:
            for line in ref_ind:
                parts = line.strip().split()
                # the 1st column is contig length
                ref_dict[parts[0]] = int(parts[1])

        # define the outfile 
        outfile = root / dirs[0] / "run.{}.sh".format(wildcards.chr)
        # win size is the size of the bin 
        block_size = params.win_size
        # define the number of bins needed fully encompass the whole chromosome
        count_to = int(ref_dict[wildcards.chr])//block_size + 1
        with open(outfile, "w") as outf:
            
            # for each bin we count through define the bin start and end
            for i in range(1, count_to):
                if i == 1:
                    beg = 1
                    end = block_size
                elif i < count_to:
                    beg = (i - 1) * block_size - 499999
                    end = i * block_size
                elif i == count_to:
                    beg = (i - 1) * block_size - 499999
                    end = int(ref_dict[wildcards.chr])

                # output the longhap commands for each bin to the shell script
                print(f"perl {longhap_ex} {input.bam} {input.vcf} {wildcards.chr}:{beg}-{end} ",
                      f"{params.max_len} {params.min_bc} {params.min_link} ",
                      f"{root}/{dirs[1]}/{wildcards.chr}_{beg}_{end}.log 2> ",
                      f"{root}/{dirs[1]}/{wildcards.chr}_{beg}_{end}.error",
                      file=outf)

            # after all the bins are taken care of
            # output the commands to merge the outputs
            # also generate het and hom and unknown barcode files
            # Finally the VCF
            print(f"perl {longhap_merge} {root}/{dirs[1]} {wildcards.chr} ",
                  f"{params.win_size} {root}/{dirs[2]}/{wildcards.chr}.log",
                  file=outf)
            print(f"cat {root/dirs[1]}/{wildcards.chr}_*.hete.barcodes | sort ",
                  f"| uniq > {root}/{dirs[2]}/{wildcards.chr}.hete.barcodes",
                  file=outf)
            print(f"cat {root/dirs[1]}/{wildcards.chr}_*.homo.barcodes | sort ",
                  f"| uniq > {root}/{dirs[2]}/{wildcards.chr}.homo.barcodes",
                  file=outf)
            print(f"cat {root/dirs[1]}/{wildcards.chr}_*.unknown.barcodes | sort",
                  f"| uniq >{root}/{dirs[2]}/{wildcards.chr}.unknown.barcodes",
                  file=outf)
            print(f"perl {longhap_log2vcf} {root}/{dirs[2]}/{wildcards.chr}.log ",
                  f"{root}/{dirs[2]}/{wildcards.chr}.hete.barcodes ",
                  f"{root}/{dirs[2]}/{wildcards.chr}.homo.barcodes ",
                  f"{input.vcf} > {root}/{dirs[2]}/{wildcards.chr}.vcf", file=outf)


# run all of the shell scripts previously defined
rule run_longhap_shells:
    input:
        "Make_Vcf/step4_longhap/LibShDir/run.{chr}.sh"
    output:
        "Make_Vcf/step4_longhap/LibOutDir/{chr}.vcf"
    shell:
        "bash {input}"


# run the longhap evaluation scripts
# the validate and cal_ori_check scripts calculate N50 and NA50 numbers
# the evaluate.pl script calculates switch errors
# it doesn't need to be run unless this is a benchmark sample
rule evaluate_longhap:
    input:
        vcfs = expand("Make_Vcf/step4_longhap/LibOutDir/{chr}.vcf", chr=CHROMS),
        vcf = "Make_Vcf/step4_longhap/{}_sentieon_pass_vars.vcf".format(config['samples']['id'])
    output:
        "Make_Vcf/step4_longhap/longhap_results.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        eval_vcfs = config['benchmark']['truth_vcf_dir'] # might be able to use the hapcut files.
    run:
        from pathlib import Path
        valid_pl = Path(params.toolsdir) / "tools" / "validate.pl"
        cal_ori = Path(params.toolsdir) / "tools" / "cal_ori_check.pl"
        eval_pl = Path(params.toolsdir) / "tools" / "evaluate.pl"

        # run the validate and cal_ori_check scripts to get N50 and NA50
        commands = [['perl', str(valid_pl), input.vcf, 'Make_Vcf/step4_longhap/LibOutDir',
                     '|', 'tee', output[0]],
                    ['perl', str(cal_ori), input.vcf, 'Make_Vcf/step4_longhap/LibOutDir',
                     '|', 'tee', '-a', output[0]]]

        # if this is a benchmark sample also run the evaluate script
        # dependant on user supplying the appropriate eval_vcfs parameter
        if params.eval_vcfs:
            commands.append(['perl', str(eval_pl), params.eval_vcfs,
                             'Make_Vcf/step4_longhap/LibOutDir', '|',
                             'tee', '-a', output[0]])

        for command in commands:
            shell(' '.join(command))

