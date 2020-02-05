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


#def determine_barcode_list(n_samples):
#    import gzip
#    import sys
#    
#    barcode_file = config['params']['toolsdir'] + config['params']['barcode']
#    barcode_rc_file = config['params']['toolsdir'] + config['params']['barcode_RC']
#    fastq_path = 
#    
#    def get_barcodes(barcodes_file):
#        print(f"reading in {barcodes_file}", file=sys.stderr)
#        nucs = ['A', 'C', 'G', 'T']
#        barcodes = []
#        with open(barcodes_file, "r") as bcs_file:
#            for line in bcs_file:
#                bc = line.strip().split()[0]
#                for i in range(0, len(bc)):
#                    for nuc in nucs:
#                        bc_alt = list(bc)   
#                        bc_alt[i] = nuc
#                        barcodes.append("".join(bc_alt))
#
#            barcodes_set = set(barcodes)
#            return barcodes_set
#
#
#    def get_fq_barcodes(fastq_path, n_samples):
#        with gzip.open(fastq_path, "rb") as fastq:
#            fq_bcs = []
#            counter = 1
#            print("Getting fastq barcodes", file=sys.stderr)
#            for line in fastq:
#                counter += 1
#                if counter == n_samples + 1:
#                    break
#                if counter % 400000 == 0:
#                    print(f"Read in {counter/4} lines", file=sys.stderr)
#                if counter % 4 == 3:
#                    seq = line.decode('utf-8').strip()
#                    fq_bcs.append(seq[-10:])
#                    fq_bcs.append(seq[-26:-16])
#                    fq_bcs.append(seq[-42:-32])
#                        
#            return fq_bcs
#
#
#    barcodes = get_barcodes(barcode_file)
#    barcodes_rc = get_barcodes(barcode_rc_file)
#    fq_bcs = get_fq_barcodes(fastq_path, n_samples)
#
#    bcs_found = 0
#    rc_bcs_found = 0
#    for bc in fq_bcs:
#        if bc in barcodes:
#            bcs_found += 1
#        if bc in barcodes_rc:
#            rc_bcs_found += 1
#
#    print(f"Barcodes: {bcs_found}\nRC_Barcodes: {rc_bcs_found}\n",file=sys.stderr)
#    if bcs_found > rc_bcs_found:
#        return config['params']['barcode']
#    elif rc_bcs_found > bcs_found:
#        return config['params']['barcode_RC']
#    else:
#        if n_samples > 40000000:
#            print("Can't determine correct barcode file for sample", file=sys.stderr)
#            sys.exit(1)
#        else:
#            determine_barcode_list(n_samples * 2)


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
