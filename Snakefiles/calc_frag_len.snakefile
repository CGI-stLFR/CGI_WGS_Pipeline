
min_frag = config['calc_frag']['min_frag']
split_dist = config['calc_frag']['split_dist']

rule remove_duplicates:
    input:
        "Align/{id}.sort.rmdup.bam"
    output:
        "Calc_Frag_Length/step1_removedup_rm000/{id}.sort.removedup_rm000.sam"
    shell:
        "samtools view -h -F 0x400 {input} | "
        "awk -F $'\t' '($1!~/#0_0_0$/){{print}}' > {output}"


rule calc_frag_len:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Calc_Frag_Length/frag_length_distribution.pdf",
        "Calc_Frag_Length/n_read_distribution.pdf",
        "Calc_Frag_Length/frag_and_bc_summary.txt",
        "Calc_Frag_Length/frag_and_bc_dataframe.tsv",
        "Calc_Frag_Length/frags_reads_per_bc.tsv",
        "Calc_Frag_Length/frags_per_bc.pdf"
    params:
        min_frag = config['calc_frag']['min_frag'],
        split_dist = config['calc_frag']['split_dist'],
        read_len = config['params']['read_len'],
        chroms = ','.join(CHROMS),
        toolsdir = config['params']['toolsdir'],
        include_dups = config['calc_frag']['include_dups']
    threads: config['threads']['calc_frag']
    run:
        if params.include_dups:
            shell("source /home/eanderson/.virtualenvs/General3/bin/activate ; "
                  "{params.toolsdir}/tools/calc_frag_len.py "
                      "--minfrag {params.min_frag} "
                      "--splitdist {params.split_dist} "
                      "--readlen {params.read_len} "
                      "--threads {threads} "
                      "--chroms {params.chroms} "
                      "--includedups "
                      "{input}")
        else:
            shell("source /home/eanderson/.virtualenvs/General3/bin/activate ; "
                  "{params.toolsdir}/tools/calc_frag_len.py "
                      "--minfrag {params.min_frag} "
                      "--splitdist {params.split_dist} "
                      "--readlen {params.read_len} "
                      "--threads {threads} "
                      "--chroms {params.chroms} "
                      "{input}")
