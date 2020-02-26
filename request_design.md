Potentially, if this could be based off of a config file at some location so that I can update the defaults without having to bother you that would be sweet!

samples:
    lanes:
        - Should autopopulate slide_lane
        - Sometimes multiple lanes from multiple slides are one sample, so that would ideally be a possibility
    id:
        - Should just be slide_lane or slide1_lane1_slide2_lane2
    chroms:
        - Should be either chr1 through chrY (as the example) if hg38 is the ref or empty if other ref
modules:
    variant_calling:
        - boolean, default false
    benchmarking:
        - boolean, default false
    phasing:
        - boolean, default false
    mate_pair_analysis:
        - boolean, default false
    read_overlap_analysis:
        - boolean, default false
    duplicate_plot:
        - boolean, default false
params:
    ref:
        - hg38 should be an option, 16s should be an option. Anything else and they can ask me ... or I can add more options when necessary
        - hg38 location: /research/rv-02/home/qmao/DBs/hg38_fa_from_NCBI/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
        - 16s location: /research/rv-02/home/eanderson/Resources_And_DBs/16sGenomes_rev.fa
    read_len:
        - ideally autopopulates based on R1
benchmark:
    - if benchmarking == True
    - options:
        - HG001/NA12878
        - HG002/NA24385
        - HG003/NA24149
        - HG004/NA24143
        - HG005/NA24631
        - HG006/NA24694
        - HG007/NA24695
    benchmark_snp:
        - "/research/rv-02/home/eanderson/Resources_And_DBs/hg38_GIAB/" + HG00# + "/benchmark_snp.chrfix.vcf.gz"
    benchmark_indel:
        - "/research/rv-02/home/eanderson/Resources_And_DBs/hg38_GIAB/" + HG00# + "/benchmark_indel.chrfix.vcf.gz"
    benchmark_bedfile:
        - "/research/rv-02/home/eanderson/Resources_And_DBs/hg38_GIAB/" + HG00# + "/benchmark_bed.chrfix.bed"

        
