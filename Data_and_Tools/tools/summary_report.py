#!/usr/bin/env python

import sys


def get_split_read_log(log_path):
    with open(log_path, "r") as log:
        for i in range(0,4):
            line = log.readline()
            if i == 1:
                real_bc = line.strip().split()[2]
            if i == 2:
                read_pair = line.strip().split()[2]
            if i == 3:
                per_split = line.strip().split()[4][1:] + "%"
                
    return real_bc, read_pair, per_split


def get_picard_metrics(log_path):
    with open(log_path, "r") as metrics:
        for line in metrics:
            if line.startswith("Mapping rate"):
                map_rate = line.strip().split()[2]
            if line.startswith("Duplicate rate"):
                dup_rate = line.strip().split()[2]
                
    return map_rate, dup_rate
                

def get_coverage_depth(log_path):
    with open(log_path, "r") as coverage:
        for i in range(0,2):
            line = coverage.readline()
            if i == 0:
                dep = line.strip().split()[3]
            if i == 1:
                cov = line.strip().split()[1]
                
    return dep, cov


def get_insert_size(log_path):
    with open(log_path, "r") as insert_size:
        for i in range(0,3):
            line = insert_size.readline()
            if i == 2:
                median = line.strip().split()[0]
                mean = line.strip().split()[4]
    
    return median, mean


def get_barcode_counts(path):
    with open(path, "r") as barcodes:
        for i in range(0,2):
            line = barcodes.readline()
            if i == 1:
                parts = tuple(line.strip().split())
                
    return parts


def get_hapcut_results(path):
    with open(path, "r") as hapcut:
        end = False
        for line in hapcut:
            if line.strip() == "combine all chrs":
                end = True
            if end == True and line.strip().startswith("phased"):
                phased = line.split()[2]
            if end == True and line.startswith("AN50:"):
                an_fifty = line.strip().split()[1]
            if end == True and line.startswith("N50:"):
                n_fifty = line.strip().split()[1]
                
    return phased, an_fifty, n_fifty
            
                
def get_longhap_results(path):
    with open(path, "r") as longhap:
        for line in longhap:
            if line.startswith("N50:"):
                n_fifty = line.strip().split(":")[1]
            if line.startswith("AN50:"):
                an_fifty = line.strip().split(":")[1]
            if line.startswith("phased_snp:"):
                phased = line.strip().split()[0].split(":")[1]
    
    return phased, an_fifty, n_fifty


if __name__ == "__main__":
    import sys
    
    try:
        sample_id = sys.argv[1]
        read_len = int(sys.argv[2])
        min_frag = sys.argv[3]
    except:
        print("Please provide a sample id, read length, and minimum fragment length")
        sys.exit(1)
    
    
    split_read_log = "split_stat_read1.log"
    picard_metrics = "Align/picard_align_metrics.txt"
    coverage_depth = "Align/coverage_depth.txt"
    insert_size_metrics = "Align/sentieon_is_{}_metric.txt".format(sample_id)
    fragment_count = "Calc_Frag_Length/step4_reassign_barcode_mq/minfrag{}_mean_len_readcount.txt".format(min_frag)
    barcode_count = "Calc_Frag_Length/step4_reassign_barcode_mq/min{}frag_barcode_summary.txt".format(min_frag)
    longhap = "Make_Vcf/step2_longhap/longhap_results.txt"
    hapcut = "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_comparison_with_giab.txt"
    
    
    real_bc, read_pair, per_split = get_split_read_log(split_read_log)
    map_rate, dup_rate = get_picard_metrics(picard_metrics)
    dep, cov = get_coverage_depth(coverage_depth)
    median, mean = get_insert_size(insert_size_metrics)
    bc_counts, frag_per_bc, bc_read_count = get_barcode_counts(barcode_count)
    frag_count, frag_len, reads_per_frag = get_barcode_counts(fragment_count)
    phased_hap, an_fifty_hap, n_fifty_hap = get_hapcut_results(hapcut)
    phased_lh, an_fifty_lh, n_fifty_lh = get_longhap_results(longhap)

    
    print("Sample ID:\t" + sample_id)
    print("Barcode Count:\t" + real_bc)
    print("Read Pairs:\t" + read_pair)
    print("Total Bases:\t" + str(int(read_pair)*read_len*2))
    print("Barcode Split Rate:\t" + per_split)
    print("Mapping Rate:\t" + map_rate)
    print("Coverage:\t" + cov)
    print("Depth:\t" + dep)
    print("Duplicate Rate:\t" + dup_rate)
    print("Median Insert Size:\t" + median)
    print("Mean Insert Size:\t" + mean)
    print("Barcode Counts:\t" + bc_counts)
    print("Average Fragment Per Barcode:\t" + frag_per_bc)
    print("Average Barcode Read Count:\t" + bc_read_count)
    print("Fragment Counts:\t" + frag_count)
    print("Average Fragment Length:\t" + frag_len)
    print("Average Reads Per Fragment:\t" + reads_per_frag)
    print("Hapcut AN50:\t" + an_fifty_hap)
    print("Hapcut N50:\t" + n_fifty_hap)
    print("Hapcut Phased Count:\t" + phased_hap)
    print("LongHap AN50:\t" + an_fifty_lh)
    print("LongHap N50:\t" + n_fifty_lh)
    print("LongHap Phased SNPs:\t" + phased_lh)
