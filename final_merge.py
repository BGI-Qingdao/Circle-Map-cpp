#!/usr/bin/env python3
import numpy as np
import pandas as pd
import pysam as ps
import pybedtools as bt
import datetime

########################################
# thresholds
########################################
fraction = 0.99
splits=0
########################################
# thresholds
########################################
def merge_fraction(chrom1,x1,x2,chrom2,y1,y2):
    """compute overlap (reciprocal) of the interval y over interval x"""
    distance = (np.minimum(x2.values,y2.values) - np.maximum(x1.values,y1.values))
    one_overlap_two = distance/(y2.values-y1.values)
    two_overlap_one = distance/(x2.values-x1.values)
    # check if they are on the same chromosome and the amount of overlap if so
    return(pd.Series(chrom1 == chrom2) + pd.Series(two_overlap_one.clip(0)) + pd.Series(one_overlap_two.clip(0)))

########################################
# load data 
########################################
raw_datas = pd.read_csv("t1", sep='\t', header=0, compression='infer', comment='#')
raw_datas.columns = ['chrom' ,'start','end','discordants','split','score']

### merge ecc intervals
norm_fraction = (fraction*2)+1
second_merging_round = raw_datas.sort_values(['chrom', 'start', 'end']).reset_index()
final_output = second_merging_round.groupby(
    merge_fraction(second_merging_round.chrom.shift(), second_merging_round.start.shift(),
                 second_merging_round.end.shift(),second_merging_round.chrom,second_merging_round.start,second_merging_round.end).lt(norm_fraction).cumsum()).agg(
    {'chrom': 'first', 'start': 'min', 'end': 'max', 'discordants' : 'max', 'split': 'sum','score':'sum'})
unfiltered_output = bt.BedTool.from_dataframe(final_output)

# filter splits
filtered = []
for interval in unfiltered_output:
    if (int(interval[4])+int(interval[3])) >= splits:
        if int(interval[1]) != 0:
            interval[1] = int(interval[1])+1
        filtered.append(interval)
# print data
filtered_output = bt.BedTool(filtered)
filtered_output.saveas("ttt1.bed")
