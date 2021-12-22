#!/bin/bash

#export PATH=/dellfsqd2/ST_OCEAN/USER/guolidong/software/bedtools2/bin:$PATH

BEDTOOLS=/dellfsqd2/ST_OCEAN/USER/guolidong/software/bedtools2/bin/bedtools
MERGEBED=/dellfsqd2/ST_OCEAN/USER/guolidong/software/bedtools2/bin/mergeBed

if [[ $# -lt 1 || $1 == "-h" ||  $1 == "--help" ]] ; then
    echo "usage : $0 <xxx.bam> <output-folder> [clustering_dist default=500]"
    exit
fi

if [[ $# -lt 1 || ! -e $1 || ( -e $2 && ! -d $2 ) ]] ; then
    echo "usage : $0 <xxx.bam> <output-folder> [clustering_dist default=500]"
    exit 1
fi

BAM=$1
OUTFD=$2

if [[ ! -e $OUTFD ]] ; then
    mkdir $OUTFD
fi
CLUSTER_DISTANCE=500
if [[ $# -gt 2 && "$3" =~ ^[+-]?[0-9]+$ ]] ; then
    CLUSTER_DISTANCE=$3
    echo "LOG: use custom CLUSTER_DISTANCE=$CLUSTER_DISTANCE"
fi

$BEDTOOLS genomecov -bg -ibam $BAM  | sort -T $OUTFD -k 1,1 -k2,2n | $MERGEBED -d $CLUSTER_DISTANCE -c 4 -o mean | sort -r -n -k 4,4 > $OUTFD/peaks.bed
