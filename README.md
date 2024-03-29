# update:
There is a pipeline for calling eccDNA based on this software and singularity, please see [panxiaoguang/Circle-Map-Cpp_Pipeline](https://github.com/panxiaoguang/Circle-Map-Cpp_Pipeline)

# Circle-Map-cpp
Circle-Map is a great tool to detect thousands of extrachromosomal circular DNA (eccDNA). Thanks to this tool, we can easily get eccDNA from Circle-Seq data!

But Circle-map realigner has low performance because it uses Python. It also has some errors in the result, especially at the start site, which will incorrectly add 1 bp.

So we **used C++ to rewrite Realiger** and also **corrected all the errors** in its result! **Welcome to try and give me some suggestion!**

## What is Circle-Map?

see https://github.com/iprada/Circle-Map and https://github.com/BGI-Qingdao/Circle-Map for details.


## Dependences

* c++ std11 compiler environment
  * g++
  * make
* htslib
* zlib
* bzip2
* xz
* python3 with libraries 
  * pysam
  * argparse
  * numpy
* bedtools
  
## Install from conda
```
conda install -c bioconda circle-map-cpp
```
Or install from source codes:

## Directly download and run a singularityCE image

For how to build an image or other details, please see [panxiaoguang/Singularity_CircleMapPlus](https://github.com/panxiaoguang/Singularity_CircleMapPlus)

```bash
wget https://github.com/panxiaoguang/Singularity_CircleMapPlus/releases/download/0.0.5/panxiaoguang-Singularity_CircleMapPlus.latest.sif
```

**how to use?**

 *If you already have a singularity program, just need to*
 
 ```bash
 ## for readExtractor
 singularity exec panxiaoguang-Singularity_CircleMapPlus.latest.sif circle_map++ ReadExtractor -i qname_unknown_circle.bam -o circular_read_candidates.bam
 ## for realign
 singularity exec --bind /home/panxiaoguang/Genome/hg38_bwa/:/hg38_bwa panxiaoguang-Singularity_CircleMapPlus.latest.sif circle_map++ Realign -i sort_circular_read_candidates.bam -qbam qname_unknown_circle.bam -sbam  sorted_unknown_circle.bam -fasta /hg38_bwa/hg38.fa -o my_unknown_circle.bed
 ## for more paramaters
 singularity exec panxiaoguang-Singularity_CircleMapPlus.latest.sif circle_map++ --help
 ```
 Note you should use `--bind` to bind a local dir into img dir
 
 *If you don't have a singularity program, you can install it or using conda version.*

## Install from source

1. you should firstly compile htslib according to this instructions : [Install](https://github.com/samtools/htslib/blob/develop/INSTALL)

2. then add the htslib to your $ENV like this:
```bash
export C_INCLUDE_PATH="your htslib path/include":$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH="your htslib path/include":$CPLUS_INCLUDE_PATH
export LIBRARY_PATH="your htslib path/lib":$LIBRARY_PATH
export LD_LIBRARY_PATH="your htslib path/lib":$LD_LIBRARY_PATH
```
3. next, make file
```bash
make
```

4. the prior pipeline will compile `realign_cm.cpp`, `merge_result.cpp` and `edlib`


**Please take note of your RAM! we load all the IO files created by the original Circle-Map into the RAM, we will use more RAM than the Circle-Map. Actually, if you have a deep sequencing of circle-seq, you will use more than 70 GB of RAM!**

## Basic usage 

```
Usage: circle_map++ <subprogram> [options]

Commands:
   ReadExtractor   Extract circular DNA read candidates
   Realign         Realign circular DNA read candidates

```
## ReadExtractor usage

```
usage: circle_map++ ReadExtractor [options]

Extracts circular DNA read candidates

required arguments:
  -i                    Input: query name sorted bam file

optional arguments:
  -o , --output         Output: Reads indicating circular DNA structural variants
  -t , --threads        Number of threads to use. Default 1
  -dir , --directory    Working directory, default is current working directory
  -q , --quality        bwa-mem mapping quality cutoff. Default value 10
  -nd, --nodiscordant   Turn off discordant (R2F1 oriented) read extraction
  -nsc, --nosoftclipped
                        Turn off soft-clipped read extraction
  -nhc, --nohardclipped
                        Turn off hard-clipped read extraction
  -v , --verbose        Verbose level, 1=error,2=warning, 3=message

```

## Realign usage

```
usage: circle_map++ Realign [options]

Realign circular DNA read candidates

Input/Output options:
  -i                    Input: bam file containing  reads extracted by ReadExtractor
  -qbam                 Input: query name sorted bam file
  -sbam                 Input: coordinate sorted bam file
  -fasta                Input: Reference genome fasta file
  -o , --output         Output filename

Running options:
  -t , --threads        Number of threads to use. Default 1
  -dir , --directory    Working directory, default will create a tmp_${pid} folder in the current working directory and automatically delete it when exit.
  -N, --no_coverage     Don't compute coverage statistics

Candidate intervals:
  -K , --clustering_dist
                        Cluster reads that are K nucleotides apart in the same node. Default: 500

Insert size estimation options:
  -ss , --sample_size   Number of concordant reads (R2F1) to use for estimating the insert size distribution. Default 100000
  -iq , --insert_mapq   Mapq cutoff for estimating the insert size distribution. Default 60

Interval processing options:
  -m , --mean_is_size   mean value of insert size
  -di , --sd_val_insert
                        SD Value of insert size
  -S , --std_factor     std_factor, extern realign interval by mIS+sIS*std_factor.(default 4)
  -q , --mapping_qual   minimum mapping quality (default 20)
  -p , --min_interval_prob
                        minimum interval probability (default 0.01)
  -e , --edit_dist_fraction
                        edit distance fraction (default 0.05)
  -l , --min_softclip_len
                        minimum softclip length (default 8)
  -rn , --max_aln_num   nhit, maximum alignment number(default 10)
  -G , --penity_gap_open
                        penalty for gap open (default 5)
  -E , --penity_gap_extern
                        penalty for gap extern (default 1)
  -P , --aln_prob       alignment probability (default 0.99)

Merge result options:
  -f , --merge_fraction
                        Merge intervals reciprocally overlapping by a fraction. Default 0.99
  -af , --allele_frequency
                        Minimum allele frequency required to report the circle interval. Default (0.1)
  -O , --number_of_discordants
                        Number of required discordant reads for intervals with only discordants. Default: 3
  -T , --split          Number of required split reads to output an eccDNA. Default: 0
  -Q , --split_quality
                        Minimum split score to output an interval. Default (0.0)
  -bs , --bases         Number of bases to extend for computing the coverage ratio. Default: 200
  -ce , --extension     Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100
  -r , --ratio          Minimum in/out required coverage ratio. Default: 0.0

```

