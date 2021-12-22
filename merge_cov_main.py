#!/usr/bin/env python3
import numpy as np
import pandas as pd
import pysam as ps
import pybedtools as bt
import datetime
import argparse
from Coverage import *
## thresholds

def merge_fraction(chrom1,x1,x2,chrom2,y1,y2):
    """compute overlap (reciprocal) of the interval y over interval x"""
    distance = (np.minimum(x2.values,y2.values) - np.maximum(x1.values,y1.values))
    one_overlap_two = distance/(y2.values-y1.values)
    two_overlap_one = distance/(x2.values-x1.values)
    # check if they are on the same chromosome and the amount of overlap if so
    return(pd.Series(chrom1 == chrom2) + pd.Series(two_overlap_one.clip(0)) + pd.Series(one_overlap_two.clip(0)))

def merge_result(interval_result1, sorted_bam, score, af, splits, n_discordant, fraction):
    ## load raw ecc intervals
    raw_datas = pd.read_csv(interval_result1, sep='\t', header=0, compression='infer', comment='#')
    raw_datas.columns = ['chrom' ,'start','end','discordant','split','score']
    raw_bed = bt.BedTool.from_dataframe(raw_datas)

    ## open soted bam of all reads
    bam = ps.AlignmentFile(sorted_bam, "rb")

    ## filter raw intervals by thresholds
    write = []
    for interval in raw_bed:
        if int(interval[4]) != 0:
            if (int(interval[4])) >= splits and float(interval[5]) > score:
                start_cov = bam.count(contig=interval[0],
                                  start=int(interval[1]), stop=int(interval[1])+1
                                  ,read_callback='nofilter')

                end_cov = bam.count(contig=interval[0],
                                start=int(interval[2])-1, stop=int(interval[2])
                                ,read_callback='nofilter')

                circle_af = ((int(interval[4]) * 2)) / ((start_cov+end_cov+0.01)/2)
                if circle_af >=af:
                    write.append(interval)
        else:
            if int(interval[3]) >= n_discordant:
                start_cov = bam.count(contig=interval[0],start=int(interval[1]), stop=int(interval[1]) + 1,
                                 read_callback='nofilter')

                end_cov = bam.count(contig=interval[0],
                                 start=int(interval[2]) - 1, stop=int(interval[2]),
                                 read_callback='nofilter')
                circle_af = (int(interval[3])) / ((start_cov+end_cov+0.01)/2)
                if circle_af >= af:
                    write.append(interval)

    bam.close()
#    for x in write:
#        print(x)
    ### merge ecc intervals
    norm_fraction = (fraction*2)+1
    unparsed_bed = bt.BedTool(write)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S:"),"Writting final output to disk")
    unparsed_pd = unparsed_bed.to_dataframe(
        names=['chrom', 'start', 'end', 'discordants', 'sc','score'])
    second_merging_round = unparsed_pd.sort_values(['chrom', 'start', 'end']).reset_index()
    #merge the output
    # merge_fraction calculates the degree of overlap between the two genomic intervals
    #lt(norm_freaction) looks the ones that surpass the merging threshold (returns 0 if true, 1 if not)
    # Cumsum calculates the cumulative sum over the output of lt. Which is then used for the grouping.
    #If the cumulative sum is the same for two rows, they are merged
    final_output = second_merging_round.groupby(
        merge_fraction(second_merging_round.chrom.shift(), second_merging_round.start.shift(),
                 second_merging_round.end.shift(),second_merging_round.chrom,second_merging_round.start,second_merging_round.end).lt(norm_fraction).cumsum()).agg(
       {'chrom': 'first', 'start': 'min', 'end': 'max', 'discordants' : 'max', 'sc': 'sum','score':'sum'})
    unfiltered_output = bt.BedTool.from_dataframe(final_output)
    # filter splits
    filtered = []
    for interval in unfiltered_output:
        if (int(interval[4])+int(interval[3])) >= splits:
            if int(interval[1]) != 0:
                interval[1] = int(interval[1])+1
            filtered.append(interval)
    filtered_output = bt.BedTool(filtered)
    return filtered_output



def coverage_caller(sorted_bam, fbed, bases, cmapq, extension, directory, ratio):
    coverage_object = coverage(sorted_bam, fbed, bases, cmapq, extension, directory)

    #Generator function for the coverage calculations
    output = coverage_object.compute_coverage(coverage_object.get_wg_coverage())
    coverage_output = filter_by_ratio(output, ratio)
    return coverage_output

'''
interval_result1 = "ttt"
sorted_bam = "./removeMT_preData/sorted_Exo_clean_hela_circle.bam"

splits = 0
score = 0.0
af = 0.1
n_discordant = 3
fraction = 0.99
directory = './'

sorted_bam = './removeMT_preData/sorted_Exo_clean_hela_circle.bam'
bases = 200
cmapq = 0
extension = 100
directory = './'
ratio= 0.0
'''

def intervalMergeCoverage(args):
    nocoverage = args.NC
    if no_coverage:
       filtered_output = merge_result(args.itab, args.sbam, args.score, args.af, args.splits, args.n_discordant, args.fraction)
       filtered_output.saveas(directory + "Exo_clean_hela_circle_site_nocov.bed")
    else:
       filtered_output = merge_result(args.itab, args.sbam, args.score, args.af, args.splits, args.n_discordant, args.fraction)
       coverage_output = coverage_caller(args.sbam, filtered_output, args.bases, args.cmapq, args.extension, args.directory, args.ratio)
       coverage_output.saveas(directory + "Exo_clean_hela_circle_site_wscov.bed")

if __name__ == '__main__':

        parser = argparse.ArgumentParser()

        parser.add_argument('-itab', '--input_table',  required=True,  type=str,    help='Input interval discordant table file')
        parser.add_argument('-sbam', '--sorted_bam',   required=True,  type=str,    help='Input sorted bam file')
##
        parser.add_argument('-dir', '--directory', metavar='',
                                 help="Working directory, default is the working directory",
                                 default=os.getcwd())
##for merge parser
        parser.add_argument('-NC', '--no_coverage', help="Don't compute coverage statistics",
                                          action='store_true')

        parser.add_argument('-S', '--split', type=int, metavar='',
                                      help="Number of required split reads to output a eccDNA. Default: 0",
                                      default=0)
        parser.add_argument('-Q', '--split_quality', type=float, metavar='',
                                           help="Minium split score to output an interval. Default (0.0)",
                                           default=0.0)
        parser.add_argument('-F', '--allele_frequency', type=float, metavar='',
                                  help="Minimum allele frequency required to report the circle interval. Default (0.1)",
                                  default=0.1)
        parser.add_argument('-O', '--number_of_discordants', type=int, metavar='',
                                      help="Number of required discordant reads for intervals with only discordants. Default: 3",
                                      default=3)
        parser.add_argument('-f', '--merge_fraction', type=float, metavar='',
                                  help="Merge intervals reciprocally overlapping by a fraction. Default 0.99",
                                  default=0.99)
##for coverage parser
        parser.add_argument('-b', '--bases', type=int, metavar='',
                                          help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                          default=200)
        parser.add_argument('-cq', '--cmapq', type=int, metavar='',
                                          help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                          default=0)
        parser.add_argument('-E', '--extension', type=int, metavar='',
                                          help="Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100",
                                          default=100)
        parser.add_argument('-r', '--ratio', type=float, metavar='',
                                      help="Minimum in/out required coverage ratio. Default: 0.0",
                                      default=0.0)

        args = vars(parser.parse_args())
        intervalMergeCoverage(args)
