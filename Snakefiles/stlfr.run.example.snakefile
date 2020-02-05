configfile: "config.yaml"


REF = config['params']['ref_fa']
if not config['samples']['chroms']:
    CHROMS = []
    try:
        with open(REF+".fai", "r") as fasta_index:
            for line in fasta_index:
                CHROMS.append(line.strip().split()[0])

    except Exception as e:
        print(f"{REF}.fai cannot be opened for parsing.")
        print(f"Exception: {e}")
        sys.exit(1)
else:
    CHROMS = config['samples']['chroms']


rule run_all:
    input:
        "Align/{}_dedup_metrics2.txt".format(config['samples']['id']),
        "Align/sentieon_metrics_{}.pdf".format(config['samples']['id']),
        "Align/flagstat_metric.txt",
        "Align/picard_align_metrics.txt",
        "Align/Distance_Between_Reads/distance_between_reads.pdf",
        "Align/Duplicate_Analysis/duplicate.pdf",
        "Align/Distance_Between_Reads/mate_pairs_distance.pdf",
        "Align/coverage_depth.txt",
        "Calc_Frag_Length/frag_and_bc_stats.txt",
        "Make_Vcf/step1_haplotyper/{}_sentieon.vcf".format(config['samples']['id']),
        "Make_Vcf/step2_benchmarking/snp_compare/summary.txt",
        "Make_Vcf/step2_benchmarking/indel_compare/summary.txt",
        "Align/Normalized_GC_Coverage.png",
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/hapblock_{}_chr1".format(config['samples']['id']),
        "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_comparison_with_giab.txt",
        "Make_Vcf/step4_longhap/longhap_results.txt",
        "summary_report.txt"


include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/calc_frag_len.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/remove_sam_analysis.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/make_vcf.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/metrics.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/stlfr.main.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/splitreads.snakefile"

shell.prefix('source /research/rv-02/home/eanderson/CGI_WGS_Pipeline/Data_and_Tools/bash_profile; ')


