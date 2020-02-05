rule cat_read_one:
    input:
        expand("{fq}/{lanes}/{lanes}_read_1.fq.gz", fq=config['samples']['fq_path'], lanes=config['samples']['lanes'])
    output:
        "data/read_1.fq.gz"
    run:
        if len(config['samples']['lanes']) > 1:
            shell("cat {input} > {output}")
        else:
            shell("ln -s ../{input} {output}")


rule cat_read_two:
    input:
        expand("{fq}/{lanes}/{lanes}_read_2.fq.gz", fq=config['samples']['fq_path'], lanes=config['samples']['lanes'])
    output:
        "data/read_2.fq.gz"
    run:
        if len(config['samples']['lanes']) > 1:
            shell("cat {input} > {output}")
        else:
            shell("ln -s ../{input} {output}")

rule split_reads:
    input:
        expand("data/read_{i}.fq.gz", i=range(1,3))
    output:
        expand("data/split_read.{i}.fq.gz", i=range(1,3))
    params:
        len = config['params']['read_len'],
        toolsdir = config['params']['toolsdir'],
        barcode = config['params']['barcode']
    run:
        shell("""
              perl {params.toolsdir}/tools/split_barcode_PEXXX_42_reads.pl \
              {params.toolsdir}{params.barcode} \
              {input} {params.len} data/split_read \
              2> data/split_stat_read.err
              """)
