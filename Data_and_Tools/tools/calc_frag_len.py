#!/usr/bin/env python3
"""calc_frag_len.py

Usage:
  calc_frag_len.py [--threads <int>] [--splitdist <int>] 
                   [--minfrag <int>] [--minreads <int>] 
                   [--readlen <int>] [--chroms <comma,separated,chroms>] 
                   [--includedups] <bampath>

Options:
  -h --help                    Show this screen.
  -m <int> --minfrag <int>     Minimum fragment size [default: 750].
  -s <int> --splitdist <int>   Distance between reads to call a new fragment [default: 300000].
  --minreads <int>             Minimum number of reads to keep a fragment [default: 4].
  --readlen <int>              Read length of dataset [default: 100].
  -n <int> --threads <int>     Number of threads to run in parallel [default: 1].
  --chroms                     Comma seperated list of chromosomes to use [default: all].
  --includedups                Include reads marked as duplicates.
  --outdir                     Specify output directory [default ./Calc_Frag_Length]
"""

# from random import sample, seed
from statistics import mean, stdev, median
import matplotlib.pyplot as plt
import multiprocessing as mp
import pandas as pd
import seaborn as sns
import pysam, argparse, os, sys
import ray

# Also SHOULD validate the arguments passed

def main():
    args = get_arguments()
    bam_path = args.bampath
    min_frag = args.minfrag
    split_dist = args.splitdist
    min_reads = args.minreads
    read_len = args.readlen
    include_dups = args.includedups
    n_threads = args.threads
    chroms = args.chroms
    dirname = args.outdir
    
    print(f"Calculating Fragment Lengths for Chroms {chroms}", file=sys.stderr)
    ray.shutdown()
    ray.init(num_cpus=n_threads, object_store_memory=100*1024*1024*1024)
    print(ray.available_resources(), file=sys.stderr)
    barcode_collection = get_reads(bam_path, split_dist, include_dups, 
                                   chroms, read_len)
    reads_per_bc_bins = write_out_barcode_summary(barcode_collection, dirname)
    barcode_collection = manipulate_df_gaps(barcode_collection, min_frag, min_reads)
    write_out_tsv_and_summary(barcode_collection, dirname, reads_per_bc_bins)
    
    

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("bampath", type=str, help="Path to bamfile for fragment length calculations.")
    parser.add_argument("-m", "--minfrag", type=int, default=750,
                        help="Minimum fragment size [default: 750].")
    parser.add_argument("-s", "--splitdist", type=int, default=300000,
                        help="Distance between reads to call a new fragment [default: 300000].")
    parser.add_argument("--minreads", type=int, default=4,
                        help="Minimum number of reads to keep a fragment [default: 4].")
    parser.add_argument("--readlen", type=int, default=100,
                        help="Read length of dataset [default: 100].")
    parser.add_argument("-n", "--threads", type=int, default=1,
                        help="Number of threads to run in parallel [default: 1].")
    parser.add_argument("--chroms", type=str,
                        help="Comma seperated list of chromosomes to use [default: all]")
    parser.add_argument("--includedups", action="store_true",
                        help="Include reads marked as duplicates.")
    parser.add_argument("--outdir", type=str, default="Calc_Frag_Length",
                        help="Specify output directory [default ./Calc_Frag_Length]")
    args = parser.parse_args()
    
    valid_chroms = get_chroms(args.bampath)
    
    if args.chroms:
        chroms = args.chroms.split(",") # validate these real quick...
        chrom_test = all(chrom in valid_chroms for chrom in chroms)
        if not chrom_test:
            parser.error("""Chroms supplied: {}
                            Chromosomes supplied must be valid chromosomes.
                            Valid Chromosomes: {}""".format(chroms, valid_chroms))
        args.chroms = chroms
    else:
        args.chroms = valid_chroms
    
    return args
    

def get_chroms(bam_path):
    chroms = []
    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
        for i in range(0, bamfile.nreferences):
            chroms.append(bamfile.get_reference_name(i))
    
    finally:
        bamfile.close()
    
    return chroms


def get_reads(bam_path, split_dist, include_dups, 
              chroms, read_len):
    tasks_list = []
    for chrom in chroms:
        object_id = p_processing_chroms.remote(bam_path, include_dups, split_dist, 
                                               chrom, read_len)
        tasks_list.append(object_id)
        
    barcode_collections = ray.get(tasks_list)
    assert len(barcode_collections) == len(chroms), f"length of done_tasks is not equal to chroms {len(done_tasks)} {len(chroms)}"
    
    barcode_collection = pd.concat(barcode_collections, ignore_index=True)
    
    return barcode_collection


@ray.remote
def p_processing_chroms(bam_path, include_dups, split_dist, 
                        chrom, read_len):
    barcode_coll = {} # Collection of barcodes to be turned into a df
    barcode_subs = {} # Collection of barcodes tracking sub-fragments
    
    def check_bc_tag(bc, chrom):
        bc_identifier = bc + "_" + chrom
        if bc_identifier not in barcode_subs:
            barcode_subs[bc_identifier] = 0

        else:
            last_pos = barcode_coll[bc_identifier + "_" + str(barcode_subs[bc_identifier])][2][-1]
            if read.reference_start - last_pos >= split_dist:
                barcode_subs[bc_identifier] += 1

        bc_tag = bc_identifier + "_" + str(barcode_subs[bc_identifier])
        return bc_tag
    
    
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    print(f"Processing Chrom: {chrom}", file=sys.stderr)
    
    for read in bamfile.fetch(chrom):
        bc = read.query_name.split("#")[1]
        if include_dups:
            read_flag = (read.flag & 0x900 == 0)
        else:
            read_flag = (read.flag & 0x1300 == 0)
        if bc != "0_0_0" and read.mapping_quality > 30 and read_flag:
            chrom = bamfile.get_reference_name(read.reference_id)
            bc_tag = check_bc_tag(bc, chrom)
            if bc_tag not in barcode_coll:
                barcode_coll[bc_tag] = [bc, chrom, [read.reference_start]]
            else:
                barcode_coll[bc_tag][2].append(read.reference_start)
    
    barcode_df = manipulate_df_prelim(barcode_coll, read_len)
    return barcode_df


def manipulate_df_gaps(test_barcodes, min_frag, min_reads):
    
    def get_read_gap(read_positions):
        gap_vals = []
        for i in range(1, len(read_positions)):
            gap_vals.append(read_positions[i] - read_positions[i-1])

        return gap_vals


    test_barcodes = test_barcodes[test_barcodes['N_Reads'] >= min_reads]
    test_barcodes = test_barcodes[test_barcodes['Frag_Length'] >= min_frag]
    test_barcodes['Frag_Gaps'] = test_barcodes['Positions'].apply(get_read_gap)
    test_barcodes['Min_Frag_Gap'] = test_barcodes['Frag_Gaps'].apply(min)
    test_barcodes['Max_Frag_Gap'] = test_barcodes['Frag_Gaps'].apply(max)
    test_barcodes['Mean_Frag_Gap'] = test_barcodes['Frag_Gaps'].apply(mean)
    return test_barcodes


def manipulate_df_prelim(barcode_collection, read_len):
    test_barcodes = pd.DataFrame.from_dict(barcode_collection, orient="index", columns=['Barcode', 'Chrom', 'Positions'])
    test_barcodes['Positions'] = test_barcodes['Positions'].apply(set).apply(list).apply(sorted)
    test_barcodes['N_Reads'] = test_barcodes['Positions'].apply(len)
    test_barcodes['Min_Pos'] = test_barcodes['Positions'].apply(min)
    test_barcodes['Max_Pos'] = test_barcodes['Positions'].apply(max)
    test_barcodes['Frag_Length'] = test_barcodes['Max_Pos'] - test_barcodes['Min_Pos'] + read_len
    return test_barcodes

    
def write_out_tsv_and_summary(test_barcodes, dirname, reads_per_bc_bins):
    print(f"Writing out summary files to {dirname}", file=sys.stderr)
    barcode_count  = test_barcodes['Barcode'].nunique()
    fragment_count = test_barcodes.shape[0]
    frags_per_bc = fragment_count/barcode_count
    avg_frag_len = test_barcodes['Frag_Length'].sum()/fragment_count
    med_frag_len = test_barcodes['Frag_Length'].median()
    avg_frag_read_count = test_barcodes['N_Reads'].sum()/fragment_count
    avg_bc_read_count = test_barcodes['N_Reads'].sum()/barcode_count
    reads_one = reads_per_bc_bins[0]
    reads_two = reads_per_bc_bins[1]
    reads_three = reads_per_bc_bins[2]
    reads_four = reads_per_bc_bins[3]
    reads_five_nine = reads_per_bc_bins[4]
    reads_ten_twenties = reads_per_bc_bins[5]
    reads_twenties_fifty = reads_per_bc_bins[6]
    reads_fifty_hundo = reads_per_bc_bins[7]
    reads_hundo_plus = reads_per_bc_bins[8]
    out_stats = dirname + "/frag_and_bc_summary.txt"
    with open(out_stats, "w") as frag_stats:
        print(f"Barcode Count:\t{barcode_count}\n"
              f"Avg Frag Per BC:\t{frags_per_bc}\n"
              f"Avg BC Read Count:\t{avg_bc_read_count}\n"
              f"Fragment Count:\t{fragment_count}\n"
              f"Avg Frag Length:\t{avg_frag_len}\n"
              f"Median Frag Length:\t{med_frag_len}\n"
              f"Avg Frag Read Count:\t{avg_frag_read_count}\n"
              f"Reads/BC (1):\t{reads_one}\n"
              f"Reads/BC (2):\t{reads_two}\n"
              f"Reads/BC (3):\t{reads_three}\n"
              f"Reads/BC (4):\t{reads_four}\n"
              f"Reads/BC (5, 10]:\t{reads_five_nine}\n"
              f"Reads/BC (10, 25]:\t{reads_ten_twenties}\n"
              f"Reads/BC (25, 50]:\t{reads_twenties_fifty}\n"
              f"Reads/BC (50, 100]:\t{reads_fifty_hundo}\n"
              f"Reads/BC (100+):\t{reads_hundo_plus}",
              file = frag_stats)
    
    test_barcodes.to_csv(dirname + "/frag_and_bc_dataframe.tsv", sep='\t', index=False)
    frag_plot = sns.distplot(test_barcodes['Frag_Length'])
    plt.savefig(dirname + "/frag_length_distribution.pdf")
    plt.clf()
    reads_plot = sns.distplot(test_barcodes['N_Reads'], kde=False)
    plt.savefig(dirname + "/n_read_distribution.pdf")
    plt.clf()
    
    
def write_out_barcode_summary(test_barcodes, dirname):
    try:
        # Create target Directory
        os.mkdir(dirname)
        print("Directory " , dirname,  " Created ") 
    except FileExistsError:
        print("Directory " , dirname,  " already exists")

    print(f"Writing out summary files to {dirname}", file=sys.stderr)
    reads_per_bc_bins=[]
    barcode_summary = (test_barcodes.groupby('Barcode')['N_Reads']
                           .agg(['sum', 'count'])
                           .reset_index()
                           .rename(columns={'sum':'Reads', 'count':'Frags'})
                       )
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 1].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 2].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 3].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 4].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 50)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 50) & (barcode_summary['Reads'] < 100)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] >= 100].Reads.sum())

    frag_per_bc = sns.distplot(barcode_summary['Frags'], kde=False, rug=True)
    plt.savefig(dirname + "/frags_per_bc.pdf")
    plt.clf()
    barcode_summary.to_csv(dirname + "/frags_reads_per_bc.tsv", sep='\t', index=False)

    return(reads_per_bc_bins)
    
    
if __name__ == "__main__":
    main()

