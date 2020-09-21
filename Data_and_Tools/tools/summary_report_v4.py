#!/usr/bin/env python3

import sys
import yaml
import re

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


def get_insert_size_gatk(log_path):
    with open(log_path, "r") as insert_size:
        for i in range(0,8):
            line = insert_size.readline()
            if i == 7:
                median = line.strip().split()[0]
                mean = line.strip().split()[5]
    
    return median, mean


def get_barcode_counts(path):
    with open(path, "r") as barcodes:
        frag_stats = []
        for line in barcodes:
            frag_stats.append(line.strip().split('\t')[1])

    return frag_stats


def get_hapcut_results(path):
    with open(path, "r") as hapcut:
        hapcut_results = {}
        for line in hapcut:
            if not line.strip():
                continue
            parts = re.split(r'\s{2,}', line.strip())
            hapcut_results[parts[0]] = parts[1]

    return hapcut_results['phased count:'], hapcut_results['AN50:'], hapcut_results['N50:']
            
                
def get_longhap_results(path):
    with open(path, "r") as longhap:
        for line in longhap:
            if line.startswith("N50:"):
                n_fifty = line.strip().split(":")[1]
                n_fifty = n_fifty.strip()
            if line.startswith("AN50:"):
                an_fifty = line.strip().split(":")[1]
                an_fifty = an_fifty.strip()
            if line.startswith("phased_snp:"):
                phased = line.strip().split()[0].split(":")[1]
    
    return phased, an_fifty, n_fifty


if __name__ == "__main__":
    import sys
    import yaml
    
    try:
        with open("config.yaml") as cnf_path:
            config = yaml.load(cnf_path, yaml.SafeLoader)
            sample_id = config['samples']['id']
            read_len = config['params']['read_len']
            split_dist = config['calc_frag']['split_dist']
    except:
        print("config.yaml not found", file=sys.stderr)

    
    split_read_log = "split_stat_read1.log"
    picard_metrics = "Align/picard_align_metrics.txt"
    coverage_depth = "Align/coverage_depth.txt"
    longhap = "Make_Vcf/step4_longhap/longhap_results.txt"
    hapcut = "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
    
    
    try:
        real_bc, read_pair, per_split = get_split_read_log(split_read_log)
        print("{:<32s}{:^20s}".format("Sample ID:", sample_id))
        print("{:<32s}{:^20s}".format("Barcode Count:", real_bc))
        print("{:<32s}{:^20s}".format("Read Pairs:", read_pair))
        print("{:<32s}{:^20d}".format("Total Bases:", int(read_pair)*read_len*2))
        print("{:<32s}{:^20s}".format("Barcode Split Rate:", per_split))
    except:
        print("Couldn't get split log stats", file=sys.stderr)

    try:
        map_rate, dup_rate = get_picard_metrics(picard_metrics)
        dep, cov = get_coverage_depth(coverage_depth)
        print("{:<32s}{:^20s}".format("Mapping Rate:", map_rate))
        print("{:<32s}{:^20s}".format("Coverage:", cov))
        print("{:<32s}{:^20s}".format("Depth:", dep))
        print("{:<32s}{:^20s}".format("Duplicate Rate:", dup_rate))
    except:
        print("Couldn't get mapping and depth stats", file=sys.stderr)

    try:
        insert_size_metrics = "Align/sentieon_is_{}_metric.txt".format(sample_id)
        median, mean = get_insert_size(insert_size_metrics)
        print("{:<32s}{:^20s}".format("Median Insert Size:", median))
        print("{:<32s}{:^20s}".format("Mean Insert Size:", mean))
    except:
        print("Couldn't get insert size stats, trying gatk file path", file=sys.stderr)
        try:
            insert_size_metrics = "Align/gatk_metrics_data.insert_size_metrics"
            median, mean = get_insert_size_gatk(insert_size_metrics)
            print("{:<32s}{:^20s}".format("Median Insert Size:", median))
            print("{:<32s}{:^20s}".format("Mean Insert Size:", mean))
        except:
            print("Couldn't get GATK insert size stats")

            

    for dist in split_dist:
        fragment_count = "Calc_Frag_Length_" + str(dist) + "/frag_and_bc_summary.txt"
        try:
            fragment_stats = get_barcode_counts(fragment_count)
            print("{:<32s}{:^20s}".format(f"{dist} Barcode Counts:", fragment_stats[0]))
            print("{:<32s}{:^20s}".format(f"{dist} Avg Frag/BC:", fragment_stats[1]))
            print("{:<32s}{:^20s}".format(f"{dist} Avg BC Read Count:", fragment_stats[2]))
            print("{:<32s}{:^20s}".format(f"{dist} Fragment Counts:", fragment_stats[3]))
            print("{:<32s}{:^20s}".format(f"{dist} Avg Frag Len:", fragment_stats[4]))
            print("{:<32s}{:^20s}".format(f"{dist} Median Frag Len:", fragment_stats[5]))
            print("{:<32s}{:^20s}".format(f"{dist} Avg Reads/Frag:", fragment_stats[6]))
        except:
            print("Couldn't get fragment stats", file=sys.stderr)

    try: 
        print("{:<32s}{:^20s}".format(f"Reads/BC (1):", fragment_stats[7]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (2):", fragment_stats[8]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (3):", fragment_stats[9]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (4):", fragment_stats[10]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (5, 10]:", fragment_stats[11]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (10, 25]:", fragment_stats[12]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (25, 50]:", fragment_stats[13]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (50, 100]:", fragment_stats[14]))
        print("{:<32s}{:^20s}".format(f"Reads/BC (100+):", fragment_stats[15]))
        print("{:<32s}{:^20s}".format(f"Total BC (1):", fragment_stats[16]))
        print("{:<32s}{:^20s}".format(f"Total BC (2):", fragment_stats[17]))
        print("{:<32s}{:^20s}".format(f"Total BC (3):", fragment_stats[18]))
        print("{:<32s}{:^20s}".format(f"Total BC (4):", fragment_stats[19]))
        print("{:<32s}{:^20s}".format(f"Total BC (5, 10]:", fragment_stats[20]))
        print("{:<32s}{:^20s}".format(f"Total BC (10, 25]:", fragment_stats[21]))
        print("{:<32s}{:^20s}".format(f"Total BC (25, 50]:", fragment_stats[22]))
        print("{:<32s}{:^20s}".format(f"Total BC (50, 100]:", fragment_stats[23]))
        print("{:<32s}{:^20s}".format(f"Total BC (100+):", fragment_stats[24]))
    except:
        print("Couldn't get fragment stats", file=sys.stderr)

    try:
        phased_hap, an_fifty_hap, n_fifty_hap = get_hapcut_results(hapcut)
        print("{:<32s}{:^20s}".format("Hapcut AN50:", an_fifty_hap))
        print("{:<32s}{:^20s}".format("Hapcut N50:", n_fifty_hap))
        print("{:<32s}{:^20s}".format("Hapcut Phased Count:", phased_hap))
    except:
        print("Couldn't get hapcut results", file=sys.stderr)
    
    try:
        phased_lh, an_fifty_lh, n_fifty_lh = get_longhap_results(longhap)
        print("{:<32s}{:^18s}".format("LongHap AN50:", an_fifty_lh))
        print("{:<32s}{:^18s}".format("LongHap N50:", n_fifty_lh))
        print("{:<32s}{:^20s}".format("LongHap Phased SNPs:", phased_lh))
    except:
        print("Couldn't get longhap results", file=sys.stderr)
