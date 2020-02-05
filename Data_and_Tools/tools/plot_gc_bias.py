'''
Author: Wenlan Tian
Date: 08/11/2016
Purpose: histogram of average intensity
'''

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import argparse


def applyHistogram(data, fig, ax, title):
    print "plotting gc average coverage ..."

    ax.plot(data['GC'], data['Average_coverage_normed'], linewidth = 1, color = 'r')

    ax.grid(True, color='black', linewidth=0.5)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 4)
    ax.set_xlabel("GC", fontsize=30)
    ax.set_ylabel("Average Coverage Normed", fontsize=30)
    ax.tick_params(labelsize=25)
    ax.set_title(title, fontsize=40)


def get_coverage_count(coverage_count_file):
    f = open(coverage_count_file, 'r')

    header = f.readline()
    info = f.readline().rstrip().split(',')
    coverage_count = float(info[2])

    return coverage_count

def plot_gc_average_coverage(gc_average_coverage_file, coverage_count_file, output_path, dataname):
    df = pd.read_csv(gc_average_coverage_file, sep='\t')

    print df

    coverage_count = get_coverage_count(coverage_count_file)

    df['Average_coverage_normed'] = df['Average_coverage'] / coverage_count

    print df

    ########################################################################################
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
    fig.suptitle("gc average coverage", fontsize=50)
    fig.tight_layout(pad = 15)

    applyHistogram(df, fig, ax, dataname)

    print "Saving to file ..."

    fig.savefig(output_path + 'GC_average_coverage_'+dataname + '.png')

    plt.close(fig)


def main(arguments):
    parser = argparse.ArgumentParser()

    parser.add_argument('--gc_average_coverage_path', dest='gc_average_coverage_dp', default='', help='Input gc average coverage file name (path)')
    parser.add_argument('--coverage_count_path', dest='coverage_count_dp', default='', help='Input coverage count file name (path)')
    parser.add_argument('--output_path', dest='output_dp', default='', help='Output path')
    parser.add_argument('--slide', dest='slide_name', default='', help='Slide name')

    args = parser.parse_args(arguments[1:])

    gc_average_coverage_file = args.gc_average_coverage_dp

    coverage_count_file = args.coverage_count_dp

    output_path = args.output_dp

    dataname = args.slide_name

    plot_gc_average_coverage(gc_average_coverage_file, coverage_count_file, output_path, dataname)
    

if __name__ == "__main__":
    main(sys.argv)
