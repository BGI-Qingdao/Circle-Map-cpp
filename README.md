# Circle-Map-cpp
rewrite part of Realign module in Circle-Map by c++

## what is Circle-Map ?

see https://github.com/iprada/Circle-Map and https://github.com/BGI-Qingdao/Circle-Map for details.


## dependences

* c++ std11 compiler environment
  * g++
  * make
* htslib
* python3 with libraries 
  * pysam
  * pybedtools
  * argparse
  * numpy
  * pandas
  
## install 

```
make
```

## usage 

```
usage: Circle-Map Realign [options]

Realign circular DNA read candidates

Input/Output options:
  -i                    Input: bam file containing the reads extracted by ReadExtractor
  -qbam                 Input: query name sorted bam file
  -sbam                 Input: coordinate sorted bam file
  -fasta                Input: Reference genome fasta file
  -o , --output         Output filename

Running options:
  -t , --threads        Number of threads to use.Default 1
  -dir , --directory    Working directory, default is the working directory
  -N, --no_coverage     Don't compute coverage statistics

Candidate intervals:
  -K , --clustering_dist
                        Cluster reads that are K nucleotides appart in the same node. Default: 500

Insert size estimation options:
  -ss , --sample_size   Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000
  -iq , --insert_mapq   Mapq cutoff for stimating the insert size distribution. Default 60
  -sd , --std           Standard deviations of the insert size to extend the intervals. Default 5

Interval processing options:
  -m , --mean_is_size   mean value of insert size
  -di , --sd_val_insert
                        SD Value of insert size
  -S , --std_factor     std_factor,extern realign interval by mIS+sIS*std_factor.(default 4)
  -q , --mapping_qual   minimum mapping quality(default 20)
  -p , --min_interval_prob
                        minimum interval probability(default 0.01)
  -e , --edit_dist_fraction
                        edit distance fraction(default 0.05)
  -l , --min_softclip_len
                        minimum softclip length(default 8)
  -rn , --max_aln_num   nhit, maximum alignment number(default 10)
  -G , --penity_gap_open
                        penity for gap open(default 5)
  -E , --penity_gap_extern
                        penity for gap extern(default 1)
  -P , --aln_prob       alignment probability(default 0.99)

check_1bp_cov:
  -af , --allele_frequency
                        Minimum allele frequency required to report the circle interval. Default (0.1)
  -O , --number_of_discordants
                        Number of required discordant reads for intervals with only discordants. Default: 3
  -T , --split          Number of required split reads to output a eccDNA. Default: 0
  -Q , --split_quality
                        Minium split score to output an interval. Default (0.0)

Merge result options:
  -f , --merge_fraction
                        Merge intervals reciprocally overlapping by a fraction. Default 0.99

Coverage metrics options:
  -bs , --bases         Number of bases to extend for computing the coverage ratio. Default: 200
  -cq , --cmapq         Minimum mapping quality treshold for coverage computation. Default: 0
  -ce , --extension     Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100
  -r , --ratio          Minimum in/out required coverage ratio. Default: 0.0

```
