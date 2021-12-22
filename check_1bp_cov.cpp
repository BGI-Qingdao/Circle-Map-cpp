//*******************std lib**********
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <map>
#include <unistd.h>
//*******************std end**********

//*******************hts lib**********
#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"
//*******************hts end**********

//*******************
static int threadnum = 8;
static float af = 0.1;
static float score = 0.0;
static int nsplit = 0;
static int ndiscordant = 3;
//*******************
struct EccTemp
{
    std::string ref;
    hts_pos_t begin;
    hts_pos_t end;
    int ndiscordant;
    int nsplit;
    float score;
    bool valid;

    void init(const std::string & line)
    {
        valid = false;
        std::istringstream ist(line);
        ist>>ref>>begin>>end>>ndiscordant>>nsplit>>score; 
    }
};

std::vector<EccTemp> cache;

void load(std::string in)
{
    std::string line ;
    std::ifstream ifs(in);
    // load header 
    std::getline(ifs,line);
    // load ecctemps
    EccTemp temp;
    while( ! std::getline(ifs,line).eof() )
    {
        temp.init(line);
        cache.push_back(temp);
    }
}

/*
 * @abstract: ternimate the program if bam is not sorted by reference
 **/
void guarentee_bam_format(sam_hdr_t * h)
{
    kstring_t tmp = { 0, 0, NULL };
    if( sam_hdr_find_line_id(h,"HD","SO","coordinate",&tmp) < 0)
    {
        hts_log_error("input bam is not sort by reference coordinate. exit... ");
        exit(1);
    }
    free(tmp.s);
}

void thread_main(int index, int threadnum, const char * bam)
{
    samFile *in = 0;
    bam1_t *b = NULL;
    bam_hdr_t * header = NULL;
    int r;

    if ((in = sam_open(bam, "rb")) == 0)
    {
         hts_log_error("Failed to open sorted bam : \"%s\" for reading", bam);
         exit(1);
    }
    header = sam_hdr_read(in);
    if (! header )
    {
        hts_log_error("Failed to load header from : \"%s\"", bam);
        exit(1);
    }
    guarentee_bam_format(header);
    hts_idx_t *idx = NULL;
    if ((idx = sam_index_load(in, bam)) == 0) 
    {
        hts_log_error("Fail to load the BAM index : \"%s\"", bam);
        exit(1);
    }
    b = bam_init1();
    //       ref           pos   count
    std::map<int, std::map<int , int > > count_cache;
    for(size_t i = index*100; i<cache.size(); i+= threadnum*100 )
    {
        for( size_t j = 0 ; j < 100 ; j++ )
        {   
            if( i+j >= cache.size()) break;
            auto & ecctemp = cache.at(i+j);
            if( ecctemp.nsplit < 1 )
            {
                if( ecctemp.ndiscordant < ndiscordant ) continue;
            }
            else 
            {
                if( ecctemp.nsplit < nsplit ||  ecctemp.score <= score ) continue; 
            }
            int tid = sam_hdr_name2tid(header,ecctemp.ref.c_str());
            int count_begin = 0, count_end =0 ;
            hts_itr_t *iter;
            if ( count_cache.find(tid) == count_cache.end() ||  count_cache[tid].find(ecctemp.begin) == count_cache[tid].end() )
            {
                if ((iter = sam_itr_queryi(idx, tid, ecctemp.begin, ecctemp.begin+1)) != 0)
                {
                     while ((r = sam_itr_next(in, iter, b)) >= 0) count_begin ++;
                     hts_itr_destroy(iter);
                }
                count_cache[tid][ecctemp.begin] = count_begin;
            }
            else
                count_begin = count_cache[tid][ecctemp.begin];
            if ( count_cache.find(tid) == count_cache.end() ||  count_cache[tid].find(ecctemp.end-1) == count_cache[tid].end() )
            {
                 if ((iter = sam_itr_queryi(idx, tid, ecctemp.end-1 , ecctemp.end )) != 0)
                 {
                      while ((r = sam_itr_next(in, iter, b)) >= 0) count_end ++;
                      hts_itr_destroy(iter);
                 }
                 count_cache[tid][ecctemp.end-1] = count_end;
            }
            else 
                 count_end = count_cache[tid][ecctemp.end-1];
            float meancov = (float(count_begin)+float(count_end)+0.01)/2.0 ;
            if( ecctemp.nsplit < 1 )
            {
                float curr_af = float(ecctemp.ndiscordant) / meancov ;
                if( curr_af >= af ) ecctemp.valid = true;
            }
            else
            {
                float curr_af = float(ecctemp.nsplit*2) / meancov ;
                if( curr_af >= af ) ecctemp.valid = true;
            }
        }
    }
    if (idx) hts_idx_destroy(idx);
    if (in) sam_close(in);
    if (header) sam_hdr_destroy(header);
    if (b) bam_destroy1(b);
}

void check_valid(const char * bam ,int threadnum)
{
    std::vector<std::thread *> workers;
    for( int i = 0; i<threadnum; i++ )
        workers.emplace_back(new std::thread(
             [bam,i,threadnum](){ 
                 thread_main(i,threadnum,bam); 
             }
       ));
    
    for( int i = 0; i<threadnum; i++ )
    {
        workers[i]->join();
        delete workers[i];
        workers[i] = NULL;
    }
    return ;
}

void print_result(const std::string output)
{
    std::ofstream ofs(output);
    ofs<<"chromesome\tstart\tend\tdiscordant\tsplit\tscore\n";
    for(size_t i = 0; i<cache.size(); i++ )
    {
        const auto & ecctemp = cache.at(i);
        if( ecctemp.valid )
           ofs<<ecctemp.ref<<'\t'
                     <<ecctemp.begin<<'\t'
                     <<ecctemp.end<<'\t'
                     <<ecctemp.ndiscordant<<'\t'
                     <<ecctemp.nsplit<<'\t'
                     <<ecctemp.score<<'\n';
    }
    ofs.close();
}

void usage()
{
    std::cerr<<"\
Usage: check_1bp_cov <-b sort.bam>\n\
                     <-i input.ecctemp.txt>\n\
                     <-o output.ecctemp.txt>\n\
                     [-t thread number(default 8)]\n\
                     [-n minimum number of split reads(default 0)]\n\
                     [-a allele frequency (default 0.1)]\n\
                     [-d minimum number of discordant reads(default 3)]\n\
                     [-s minimum score (default 0.0)]"<<std::endl;
}

int main(int argc, char **argv)
{
    std::string bam;
    std::string input;
    std::string output;

    int option;
    while((option = getopt(argc, argv, "b:i:o:t:n:a:d:s:h")) != -1)
    {
        switch(option)
        {
            case 'b': bam = std::string(optarg);break;
            case 'i': input = std::string(optarg);break;
            case 'o': output = std::string(optarg);break;
            case 't': threadnum = atoi(optarg);break;
            case 'n': nsplit = atoi(optarg);break;
            case 'a': af = atof(optarg);break;
            case 'd': ndiscordant = atoi(optarg);break;
            case 's': score = atof(optarg);break;
            case 'h': usage();return 0;
            default:  usage();return 1;
        }
    }
    if( bam.empty() || input.empty() || output.empty() )
    {
        usage();
        return 1;
    }
    load(input);
    check_valid(bam.c_str(),threadnum);
    print_result(output);
    return 0;
}
