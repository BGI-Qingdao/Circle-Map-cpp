#!/usr/bin/env python3
import os,sys,time
import numpy as np
import pandas as pd
import pysam as ps
import pybedtools as bt
import datetime
import argparse
from Coverage import *

pwd_config_file = os.path.realpath(__file__)
#print(pwd_config_file)
config_tool_path = '/'.join(pwd_config_file.split('/')[:-1]) 
#print(config_tool_path)
def is_soft_clipped(read):

    """Function that checks the CIGAR string of the sam file and returns true if the read is soft-clipped"""

    # cigar 4 equals to S in pysam sam representation
    match = 0
    for cigar in read.cigar:
        if cigar[0] == 4:
            match +=1
        else:
            pass

    if match > 0:
        return(True)

    else:
        return(False)

def is_hard_clipped(read):

    """Function that checks the CIGAR string of the sam file and returns true if the read is hard-clipped"""

    # cigar 5 equals to H in pysam sam representation
    match = 0
    for cigar in read.cigar:
        if cigar[0] == 5:
            match += 1
        else:
            pass

    if match > 0:
        return (True)

    else:
        return (False)

def insert_size_dist(sample_size,mapq_cutoff,qname_bam):
    """Function that takes as input a queryname sorted bam and computes the mean insert a size and
    the standard deviation from. This number is computed from the F1R2 read with a user defined sample size,
     using a user defined mapping quality cutoff in both reads"""

    whole_bam = ps.AlignmentFile(qname_bam, "rb")

    counter = 0
    insert_length = []
    read1 = ''

    # this is similar to the code of read extractor. I save the first read in memory and then I operate
    # in both reads together
    for read in whole_bam:

        if read.is_read1:
            read1 = read
        else:
            if read.is_read2 and read.qname == read1.qname:
                read2 = read
                # both reads in memory
                if read1.mapq >= mapq_cutoff and read2.mapq >= mapq_cutoff:
                    if read1.is_proper_pair:
                        if is_hard_clipped(read1) == False and is_hard_clipped(read2) == False:
                            if is_soft_clipped(read1) == False and is_soft_clipped(read2) == False:
                                if read1.is_reverse == False and read2.is_reverse == True:
                                    if read1.tlen > 0:
                                        insert_length.append(read1.tlen)
                                        counter += 1

        if counter >= sample_size:
            break
        else:
            pass


    mean = np.mean(insert_length)
    std = np.std(insert_length)
    print(mean, std)
    return(mean, std)


def merge_fraction(chrom1,x1,x2,chrom2,y1,y2):
    """compute overlap (reciprocal) of the interval y over interval x"""
    distance = (np.minimum(x2.values,y2.values) - np.maximum(x1.values,y1.values))
    one_overlap_two = distance/(y2.values-y1.values)
    two_overlap_one = distance/(x2.values-x1.values)
    # check if they are on the same chromosome and the amount of overlap if so
    return(pd.Series(chrom1 == chrom2) + pd.Series(two_overlap_one.clip(0)) + pd.Series(one_overlap_two.clip(0)))

def merge_result(interval_result1, split, fraction):
    ## load raw ecc intervals
    raw_datas = pd.read_csv(interval_result1, sep='\t', header=0, compression='infer', comment='#')
    raw_datas.columns = ['chrom' ,'start','end','discordants','sc','score']
    ### merge ecc intervals
    norm_fraction = (fraction*2)+1
    #unparsed_bed = bt.BedTool(write)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S:"),"Writting final output to disk")
    unparsed_pd = raw_datas
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
        if (int(interval[4])+int(interval[3])) >= split:
            if int(interval[1]) != 0:
                interval[1] = int(interval[1])+1
            filtered.append(interval)
    filtered_output = bt.BedTool(filtered)
    return filtered_output


def coverage_metrics(sorted_bam, fbed, bases, cmapq, extension, directory, ratio):
    coverage_object = coverage(sorted_bam, fbed, bases, cmapq, extension, directory)

    #Generator function for the coverage calculations
    output = coverage_object.compute_coverage(coverage_object.get_wg_coverage())
    coverage_output = filter_by_ratio(output, ratio)
    return coverage_output


def candidate_intervals(args):
    """find candidate intervals"""
    os.system("sh " + config_tool_path + "/candidate_intervals.sh " + ' '.join([
        args.i, 
        args.directory, 
        str(args.clustering_dist)])
    )

def insert_size_est(args):
    mean, std=insert_size_dist(args.sample_size, args.insert_mapq, args.qbam)
    return(mean, std)

def realign_cm(args, mean, std):
    os.system(config_tool_path + "/realign_cm " + ' '.join([str(i) for i in [
        '-c',args.i,
        '-b',args.directory + "/peaks.bed",
        '-g',args.fasta,
        '-o',args.directory + "/ecctemp.txt",
        '-m',mean,
        '-d',std,
        '-S',args.std_factor,
        '-t',args.threads,
        '-q',args.mapping_qual,
        '-p',args.min_interval_prob,
        '-e',args.edit_dist_fraction,
        '-l',args.min_softclip_len,
        '-n',args.max_aln_num,
        '-G',args.penity_gap_open,
        '-E',args.penity_gap_extern,
        '-P',args.aln_prob,
        ]])
    )

def check_1bp_cov(args):
    os.system(config_tool_path + "/check_1bp_cov " + ' '.join([str(i) for i in[
        '-b', args.sbam,
        '-i', args.directory + "/ecctemp.txt",
        '-o', args.directory + "/ecctemp.csv",
        '-t', args.threads,
        '-n', args.split,
        '-a', args.allele_frequency,
        '-d', args.number_of_discordants,
        '-s', args.split_quality,
        ]])
    )

def merge_nocov_result(args):
    filtered_output = merge_result(
        args.directory + "/ecctemp.csv", 
        args.split, 
        args.merge_fraction)

    filtered_output.saveas(args.output)

def coverage_caller(args):
    filtered_output = merge_result(args.directory + "/ecctemp.csv",  args.split,  args.merge_fraction) 
    coverage_output = coverage_metrics(args.sbam, filtered_output, args.bases, args.cmapq, args.extension, args.directory, args.ratio) 
    coverage_output.to_csv(args.output, header=None, index=None, sep='\t', mode='w')


class circle_map_cpp:

    #def __getpid__(self):

    #    pid = os.getpid()
    #    return (pid)

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Circle-Map++Realign',
            usage='''Circle-Map++Realign <subprogram> [options]

version=1.0.0
contact= https://github.com/BGI-Qingdao/Circle-Map++/issues

The Circle-Map++Realign

Commands:
   Realign         Realign circular DNA read candidates

''' )

        subparsers = self.parser.add_subparsers()

        self.realigner = subparsers.add_parser(
            name="Realign",
            description='Realign circular DNA read candidates',
            prog="Circle-Map Realign",
            usage='''Circle-Map Realign [options]'''
            )

        if len(sys.argv) <= 1:
            self.parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo argument given to Circle-Map++Realign"
                                 "\nExiting\n")
            sys.exit(0)

        else:
            if sys.argv[1] == "Realign":
                self.subprogram = self.args_realigner()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                candidate_intervals(self.args)
                mean, std=insert_size_est(self.args)
                realign_cm(self.args, mean, std)
                check_1bp_cov(self.args)

                if self.args.no_coverage == False :
                    coverage_caller(self.args)
                else:
                    merge_nocov_result(self.args)
                    

            else:
                self.parser.print_help()
                time.sleep(0.01)
                sys.stderr.write("\nWrong argument given to Circle-Map"
                                     "\nExiting\n")
                sys.exit(0)

    def args_realigner(self):
        parser = self.realigner

        # declare the different groups for the parser
        parser._action_groups.pop()
        io_options = parser.add_argument_group('Input/Output options')
        running = parser.add_argument_group('Running options')
        candidate_intervals = parser.add_argument_group('Candidate intervals')
        insert_size_est = parser.add_argument_group('Insert size estimation options')
        realign_cm = parser.add_argument_group('Interval processing options')
        check_1bp_cov = parser.add_argument_group('check_1bp_cov')
        merge_result = parser.add_argument_group('Merge result options')
        coverage_caller = parser.add_argument_group('Coverage metrics options')

        io_options.add_argument('-i', metavar='', help="Input: bam file containing the reads extracted by ReadExtractor")
        io_options.add_argument('-qbam', metavar='', help="Input: query name sorted bam file")
        io_options.add_argument('-sbam', metavar='', help="Input: coordinate sorted bam file")
        io_options.add_argument('-fasta', metavar='', help="Input: Reference genome fasta file")

        if "-i" in sys.argv  and "-qbam" in sys.argv and "-fasta" in sys.argv and "-sbam" in sys.argv:
            io_options.add_argument('-o', '--output', metavar='', help="Output filename",
                                    default="circle_%s.bed" % sys.argv[sys.argv.index("-i") + 1])

        # candidate intervals
            candidate_intervals.add_argument('-K', '--clustering_dist', type=int,   default=500, metavar='',   
                                                          help="Cluster reads that are K nucleotides appart in the same node. Default: 500" )
        # check coverage
            # check_1bp_cov.add_argument('-t', '--thread',                type=int,   default=8,      metavar='', 
            #                                   help="thread number(default 8)"
            #                                   )
            check_1bp_cov.add_argument('-af', '--allele_frequency',     type=float, default=0.1,    metavar='',
                                               help="Minimum allele frequency required to report the circle interval. Default (0.1)" 
                                               )
            check_1bp_cov.add_argument('-O', '--number_of_discordants',  type=int,       default=3, metavar='',
                                               help="Number of required discordant reads for intervals with only discordants. Default: 3"
                                               )
            check_1bp_cov.add_argument('-T', '--split',                  type=int,       default=0, metavar='', 
                                               help="Number of required split reads to output a eccDNA. Default: 0" 
                                               )
            check_1bp_cov.add_argument('-Q', '--split_quality',          type=float,     default=0.0, metavar='',
                                               help="Minium split score to output an interval. Default (0.0)" 
                                               )
        # insert size
            insert_size_est.add_argument('-ss', '--sample_size',        type=int,   default=100000, metavar='',
                                                 help="Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000"
                                                 )
            insert_size_est.add_argument('-iq', '--insert_mapq',        type=int,   default=60,     metavar='',
                                                 help="Mapq cutoff for stimating the insert size distribution. Default 60"
                                                 )
            realign_cm.add_argument('-m', '--mean_is_size',          type=int,      metavar='', 
                                             help="mean value of insert size" 
                                             )
            realign_cm.add_argument('-di', '--sd_val_insert',         type=int,      metavar='', 
                                             help="SD Value of insert size" 
                                             )
            realign_cm.add_argument('-S', '--std_factor',            type=int,      default=4, metavar='', 
                                             help=" std_factor,extern realign interval by mIS+sIS*std_factor.(default 4)" 
                                             )
            # realign_cm.add_argument('-t', '--num_thread',            type=int,      default=8,      metavar='', 
            #                                  help="thread number(default 8)" 
            #                                  )
            realign_cm.add_argument('-q', '--mapping_qual',          type=float,    default=20,     metavar='', 
                                             help="minimum mapping quality(default 20)" 
                                             )
            realign_cm.add_argument('-p', '--min_interval_prob',     type=float,    default=0.01,   metavar='', 
                                             help="minimum interval probability(default 0.01)" 
                                             )
            realign_cm.add_argument('-e', '--edit_dist_fraction',    type=float,    default=0.05,   metavar='', 
                                             help="edit distance fraction(default 0.05)" 
                                             )
            realign_cm.add_argument('-l', '--min_softclip_len',      type=float,    default=8,      metavar='', 
                                             help="minimum softclip length(default 8)" 
                                             )
            realign_cm.add_argument('-rn', '--max_aln_num',           type=float,    default=10,     metavar='',  
                                             help="nhit, maximum alignment number(default 10)" 
                                             )
            realign_cm.add_argument('-G', '--penity_gap_open',       type=float,    default=5,      metavar='', 
                                             help="penity for gap open(default 5)" 
                                             )
            realign_cm.add_argument('-E', '--penity_gap_extern',     type=float,    default=1,      metavar='',  
                                             help="penity for gap extern(default 1)" 
                                             )
            realign_cm.add_argument('-P', '--aln_prob',              type=float,    default=0.99,   metavar='', 
                                             help="alignment probability(default 0.99)" 
                                             )

        # merge 

            merge_result.add_argument('-f', '--merge_fraction',         type=float,     default=0.99, metavar='',
                                               help="Merge intervals reciprocally overlapping by a fraction. Default 0.99" 
                                               )
        # coverage metrics
            coverage_caller.add_argument('-bs', '--bases',               type=int,       default=200,    metavar='',
                                                  help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                                  )
            coverage_caller.add_argument('-cq', '--cmapq',              type=int,       default=0,      metavar='',
                                                  help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                                  )
            coverage_caller.add_argument('-ce', '--extension',           type=int,       default=100,    metavar='',
                                                  help="Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100",
                                                  )
            coverage_caller.add_argument('-r', '--ratio',               type=float,     default=0.0,    metavar='',
                                                  help="Minimum in/out required coverage ratio. Default: 0.0" 
                                                  )

        # run options
            running.add_argument('-t', '--threads', type=int, default=1, metavar='',
                                         help="Number of threads to use.Default 1",
                                         )
            running.add_argument('-dir', '--directory', default=os.getcwd(), metavar='',
                                         help="Working directory, default is the working directory",
                                         )
            running.add_argument('-N', '--no_coverage', help="Don't compute coverage statistics",  
                              action='store_true')                                              

        else:
            io_options.add_argument('-o', '--output', metavar='', help="Output filename")

        # candidate intervals
            candidate_intervals.add_argument('-K', '--clustering_dist', type=int,   default=500, metavar='',   
                                                          help="Cluster reads that are K nucleotides appart in the same node. Default: 500" )
        # check coverage
            check_1bp_cov.add_argument('-af', '--allele_frequency',     type=float, default=0.1,    metavar='',
                                               help="Minimum allele frequency required to report the circle interval. Default (0.1)" 
                                               )
            check_1bp_cov.add_argument('-O', '--number_of_discordants',  type=int,       default=3, metavar='',
                                               help="Number of required discordant reads for intervals with only discordants. Default: 3"
                                               )
            check_1bp_cov.add_argument('-T', '--split',                  type=int,       default=0, metavar='', 
                                               help="Number of required split reads to output a eccDNA. Default: 0" 
                                               )
            check_1bp_cov.add_argument('-Q', '--split_quality',          type=float,     default=0.0, metavar='',
                                               help="Minium split score to output an interval. Default (0.0)" 
                                               )
        # insert size
            insert_size_est.add_argument('-ss', '--sample_size',        type=int,   default=100000, metavar='',
                                                 help="Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000"
                                                 )
            insert_size_est.add_argument('-iq', '--insert_mapq',        type=int,   default=60,     metavar='',
                                                 help="Mapq cutoff for stimating the insert size distribution. Default 60"
                                                 )
            insert_size_est.add_argument('-sd', '--std',                type=int,   default=4,      metavar='',
                                                 help="Standard deviations of the insert size to extend the intervals. Default 5" 
                                                 )
        # realignment
            # realign_cm.add_argument('-c', '--candidate.bam',         type=str,      metavar='', 
            #                                  help="candidate.bam" 
            #                                  )
            realign_cm.add_argument('-m', '--mean_is_size',         type=int,      metavar='', 
                                             help="mean value of insert size" 
                                             )
            realign_cm.add_argument('-di', '--sd_val_insert',       type=int,      metavar='', 
                                             help="SD Value of insert size" 
                                             )
            realign_cm.add_argument('-S', '--std_factor',           type=int,      metavar='', 
                                             help=" std_factor,extern realign interval by mIS+sIS*std_factor.(default 4)" 
                                             )
            # realign_cm.add_argument('-t', '--num_thread',            type=int,      default=8,      metavar='', 
            #                                  help="thread number(default 8)" 
            #                                  )
            realign_cm.add_argument('-q', '--mapping_qual',         type=float,    default=20,     metavar='', 
                                             help="minimum mapping quality(default 20)" 
                                             )
            realign_cm.add_argument('-p', '--min_interval_prob',    type=float,    default=0.01,   metavar='', 
                                             help="minimum interval probability(default 0.01)" 
                                             )
            realign_cm.add_argument('-e', '--edit_dist_fraction',   type=float,    default=0.05,   metavar='', 
                                             help="edit distance fraction(default 0.05)" 
                                             )
            realign_cm.add_argument('-l', '--min_softclip_len',     type=float,    default=8,      metavar='', 
                                             help="minimum softclip length(default 8)" 
                                             )
            realign_cm.add_argument('-rn', '--max_aln_num',         type=float,    default=10,     metavar='',  
                                             help="nhit, maximum alignment number(default 10)" 
                                             )
            realign_cm.add_argument('-G', '--penity_gap_open',      type=float,    default=5,      metavar='', 
                                             help="penity for gap open(default 5)" 
                                             )
            realign_cm.add_argument('-E', '--penity_gap_extern',    type=float,    default=1,      metavar='',  
                                             help="penity for gap extern(default 1)" 
                                             )
            realign_cm.add_argument('-P', '--aln_prob',             type=float,    default=0.99,   metavar='', 
                                             help="alignment probability(default 0.99)" 
                                             )

        # merge 

            merge_result.add_argument('-f', '--merge_fraction',         type=float,     default=0.99,       metavar='',
                                               help="Merge intervals reciprocally overlapping by a fraction. Default 0.99" 
                                               )
        # coverage metrics
            coverage_caller.add_argument('-bs', '--bases',              type=int,       default=200,        metavar='',
                                                  help="Number of bases to extend for computing the coverage ratio. Default: 200",
                                                  )
            coverage_caller.add_argument('-cq', '--cmapq',              type=int,       default=0,          metavar='',
                                                  help="Minimum mapping quality treshold for coverage computation. Default: 0",
                                                  )
            coverage_caller.add_argument('-ce', '--extension',          type=int,       default=100,        metavar='',
                                                  help="Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100",
                                                  )
            coverage_caller.add_argument('-r', '--ratio',               type=float,     default=0.0,        metavar='',
                                                  help="Minimum in/out required coverage ratio. Default: 0.0" 
                                                  )

        # run options
            running.add_argument('-t', '--threads', type=int, default=1, metavar='',
                                         help="Number of threads to use.Default 1",
                                         )
            running.add_argument('-dir', '--directory', default=os.getcwd(), metavar='',
                                         help="Working directory, default is the working directory",
                                         )
            running.add_argument('-N', '--no_coverage', help="Don't compute coverage statistics",  
                              action='store_true')                                              
           # find out which arguments are missing

            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write("\nInput does not match. Check that you provide the -i, -sbam, -qbam and -fasta options"
                             "\nExiting\n")
            sys.exit(0)
            if len(sys.argv[2:]) == 0:
                parser.print_help()
                time.sleep(0.01)
                sys.stderr.write("\nNo arguments given to Realign. Exiting\n")
                sys.exit(0)

        return (parser)

def main():
    run = circle_map_cpp()

if __name__ == '__main__':
    main()
