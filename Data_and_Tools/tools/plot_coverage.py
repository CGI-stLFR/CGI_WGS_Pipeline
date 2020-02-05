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
import os.path
from scipy.stats import poisson
import argparse


def poisson_probability(actual, mean):
    # naive:   math.exp(-mean) * mean**actual / factorial(actual)

    # iterative, to keep the components from getting too large or small:
    p = math.exp(-mean)
    for i in xrange(int(actual)):
        p *= mean
        p /= i+1
    return p


def applyHistogram_relative_poisson(data, fig, ax, title, poisson):
    print "plotting histogram poisson ..."
    ax.bar(data['coverage'], data['count_ratio'],  alpha=1, linewidth = 0.1, color = 'r')
    ax.plot(data['coverage'], data[poisson], color = 'blue')

    ax.grid(True, color='black', linewidth=0.5)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 0.12)
    ax.set_xlabel("Coverage", fontsize=30)
    ax.set_ylabel("Relative Frequency", fontsize=30)
    ax.tick_params(labelsize=25)
    ax.set_title(title, fontsize=40)


def process_data_mode(data):

    data_processed = data.split('.')[-2] + '_processed_mode.csv'
    print data_processed

    if os.path.isfile(data_processed):
        print 'exist'
        df = pd.read_csv(data_processed, sep=' ')

    else:
        print 'not exist'
        df = pd.read_csv(data, sep=' ')

        total_coverage = df['count'].sum()

        #print type(pdf)

        df['count_ratio'] = df['count']/total_coverage

        df['coverage_ratio'] = df['coverage'] * df['count_ratio']

        sum_coverage_ratio = df['coverage_ratio'].sum()

        print sum_coverage_ratio

        coverage_mode= df['count'].idxmax()

        print coverage_mode

        #df['poisson'] = df.apply(lambda x: poisson.pmf(x['coverage'], sum_coverage_ratio), axis=1)

        df['poisson_1'] = df.apply(lambda x: poisson_probability(x['coverage'], coverage_mode-1), axis=1)
        df['poisson_0.5'] = df.apply(lambda x: poisson_probability(x['coverage'], coverage_mode-0.5), axis=1)
        df['poisson'] = df.apply(lambda x: poisson_probability(x['coverage'], coverage_mode), axis=1)
        df['poisson+0.5'] = df.apply(lambda x: poisson_probability(x['coverage'], coverage_mode+0.5), axis=1)
        df['poisson+1'] = df.apply(lambda x: poisson_probability(x['coverage'], coverage_mode+1), axis=1)

        #df['kernal'] = pdf

        df.to_csv(data_processed, sep=' ', index=False)

    return df


def calculate_diff(df, col):
    diff1 = 0.0
    diff2 = 0.0

    for index, row in df.iterrows():
        if row['coverage'] > 0 and row['count_ratio'] <= row[col]:
            break
        diff1 += row['count_ratio'] - row[col]

    diff1 = round(diff1, 4)
    print diff1


    for index, row in df[::-1].iterrows():
        if row['count_ratio'] <= row[col]:
            break
        diff2 += row['count_ratio'] - row[col]

    diff2 = round(diff2, 4)
    print str(diff1) + ',' + str(diff2) + ',' + str(diff1 + diff2)

    return diff1, diff2

def calculate_diff_multiple(df):
    diff1_1, diff1_2 = calculate_diff(df, 'poisson_1')
    diff2_1, diff2_2 = calculate_diff(df, 'poisson_0.5')
    diff3_1, diff3_2 = calculate_diff(df, 'poisson')
    diff4_1, diff4_2 = calculate_diff(df, 'poisson+0.5')
    diff5_1, diff5_2 = calculate_diff(df, 'poisson+1')

    diffs = [(diff1_1, diff1_2), (diff2_1, diff2_2), (diff3_1, diff3_2), (diff4_1, diff4_2), (diff5_1, diff5_2)]

    poissons = ['poisson_1','poisson_0.5', 'poisson', 'poisson+0.5', 'poisson+1']

    diff_sums = [0] * 5

    diff_sums[0] = diff1_1 + diff1_2;
    diff_sums[1] = diff2_1 + diff2_2;
    diff_sums[2] = diff3_1 + diff3_2;
    diff_sums[3] = diff4_1 + diff4_2;
    diff_sums[4] = diff5_1 + diff5_2;

    min_idx = 0
    min_sum = diff_sums[0]

    for i in range(1,5):
        if(diff_sums[i] < min_sum):
            min_sum = diff_sums[i]
            min_idx = i

    return diffs[min_idx], poissons[min_idx] 

def calculate_low_coverage(df, coverage_thredshold):

    low_coverage = 0

    for coverage in range(coverage_thredshold):
        low_coverage += df[df['coverage']==coverage]['count_ratio'].values[0]

    low_coverage = round(low_coverage, 8)


    return low_coverage


def plot_coverage_mode(coverage_file_path, output_path, dataname):
    df = process_data_mode(coverage_file_path)

    (diff1,diff2),poisson = calculate_diff_multiple(df)

    low_coverage1 = calculate_low_coverage(df, 5)

    low_coverage3 = calculate_low_coverage(df, 10)


    ########################################################################################
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
    fig.suptitle("coverage histogram", fontsize=50)
    fig.tight_layout(pad = 15)


    max_count = df['count'].max()

    print max_count

    applyHistogram_relative_poisson(df, fig, ax, dataname, poisson)

    fig.text(0.20,0.1,'left:'+str(diff1), fontsize=20)
    fig.text(0.45,0.1,'right:' + str(diff2), fontsize=20)
    fig.text(0.20,0.05,'<5x:'+str(low_coverage1), fontsize=20)
    fig.text(0.45,0.05,'<10x:'+str(low_coverage3), fontsize=20)

    #ax.grid(True, color='black', linewidth=0.5)

    print "Saving to file ..."

    fig.savefig(output_path + 'Coverage_Histogram_'+dataname + '.png')
    file_out = output_path + 'Coverage_Summary_'+dataname + '.csv'

    plt.close(fig)

    f_out = open(file_out, 'w')
    f_out.write('left,right,left+right,<5x,<10\n')
    f_out.write(str(diff1) + ',' + str(diff2) + ',' + str(diff1+diff2) + ',' + str(low_coverage1) + ',' + str(low_coverage3) + '\n')

def main(arguments):
    parser = argparse.ArgumentParser()

    parser.add_argument('--coverage_path', dest='coverage_dp', default='', help='Input coverage file name (path)')
    parser.add_argument('--output_path', dest='output_dp', default='', help='Output path')
    parser.add_argument('--slide', dest='slide_name', default='', help='Slide name')

    args = parser.parse_args(arguments[1:])

    coverage_file = args.coverage_dp

    output_path = args.output_dp

    dataname = args.slide_name

    #path = '/home/wtian/Projects/Bioinformatics/coverage_database/reference_index/output/CL100071428_PE100_L01_30x.csv'

    plot_coverage_mode(coverage_file, output_path, dataname)
    

if __name__ == "__main__":
    main(sys.argv)
