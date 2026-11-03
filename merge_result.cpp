#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#include "utils/MultiThread.h"
//*******************
static int threadnum = 8;
static float af = 0.1;
static float score = 0.0;
static int nsplit = 0;
static int ndiscordant = 3;
static float merge_factor = 0.99;
static int ext = 200;
static int ilen = 100;
static float ratio = 0.0;
//*******************

struct PreDepth
{
    //        chrname        covs
    std::map<std::string,std::vector<int> > datas;
    std::map<std::string,int> ref_lens;
    void init_fai(const std::string & genome)
    {
        std::string fai_name = genome + ".fai";
        std::ifstream ifs(fai_name);
        if( ! ifs.good() )
        {
            std::cerr<<"fai file of "<<genome<<" not exist ! exit ..."<<std::endl;
            exit(1);
        }
        std::string line;
        std::string chr_name;
        int chr_len;
        while( ! std::getline(ifs,line).eof() )
        {
            std::istringstream iss(line);
            iss>>chr_name>>chr_len;
            ref_lens[chr_name]=chr_len;
            datas[chr_name].resize(chr_len,0);      
        }
    }
    void init_covs(const std::string & fcov)
    {
        std::ifstream ifs(fcov);
        if( ! ifs.good() )
        {
            std::cerr<<"coverage not exist ! exit ..."<<std::endl;
            exit(1);
        }
        std::string line;
        std::string chr_name;
        int start;
        int end;
        int cov;
        while( ! std::getline(ifs,line).eof() )
        {
            std::istringstream iss(line);
            iss>>chr_name>>start>>end>>cov;
            auto itr = datas.find(chr_name);
            if(cov > 0 && itr != datas.end() )
            {
                for(int i = start ; i<end; i++)
                    itr->second[i] = cov;
            }
        }
    }
    int cov_at(const std::string &chr_name, int pos)
    {
        auto itr = datas.find(chr_name);
        if ( itr == datas.end() || pos <0 || pos >= ref_lens[chr_name] )
             return 0;
        return itr->second.at(pos);
    }
    std::vector<int> cov_region(const std::string &chr_name, int start, int end)
    {
        auto itr = datas.find(chr_name);
        if ( itr == datas.end())
             return std::vector<int>();
        if( start < 0 )  start = 0;
        if( end >= ref_lens[chr_name] ) end = ref_lens[chr_name]-1;
        std::vector<int> ret ;
        ret.resize(end-start+1);
        for(int i = 0 ; i < end-start+1 ;i++)
           ret[i] = itr->second[start+i];
        return ret;
    }
} cov_cache;

struct eccCache
{
    struct EccTemp
    {
        int begin;
        int end;
        int ndiscordant;
        int nsplit;
        float score;
        bool valid;
        bool operator < (const EccTemp& o)const
        {
            return begin < o.begin || ( begin == o.begin && end < o.end );
        }
        void init(const std::string & line,std::string &ref)
        {
            valid = false;
            std::istringstream ist(line);
            ist>>ref>>begin>>end>>ndiscordant>>nsplit>>score;
        }
    };

    /*
     * @abstract: true if a and b overlap reciprocally by more than
     *            merge_factor -- the same criterion the previous chained
     *            merge applied to consecutive candidates, but stated for an
     *            unordered pair.
     */
    static bool reciprocal_overlap(const EccTemp & a, const EccTemp & b)
    {
        int max_start = a.begin > b.begin ? a.begin : b.begin ;
        int min_end   = a.end   < b.end   ? a.end   : b.end   ;
        if( min_end <= max_start ) return false;
        double len_a = double(a.end - a.begin);
        double len_b = double(b.end - b.begin);
        if( len_a <= 0.0 || len_b <= 0.0 ) return false;
        double ov = double(min_end - max_start);
        return ( ov/len_a + ov/len_b ) > 2.0*merge_factor ;
    }

    /* union-find over candidate indexes */
    static size_t uf_find(std::vector<size_t> & parent, size_t x)
    {
        while( parent[x] != x )
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    static void uf_unite(std::vector<size_t> & parent, size_t a, size_t b)
    {
        size_t ra = uf_find(parent,a), rb = uf_find(parent,b);
        if( ra == rb ) return ;
        if( ra < rb ) parent[rb] = ra; else parent[ra] = rb;
    }

    std::map<std::string,std::vector<EccTemp>> temp_caches;
    
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
            std::string ref;
            temp.init(line,ref);
            temp_caches[ref].push_back(temp);
        }
        for(auto & pair :temp_caches)
        {
            auto & cache = pair.second;
            std::sort(cache.begin(),cache.end());
        }
    }
    
    void check_1bp_cov()
    {
        for(auto & pair :temp_caches)
        {
            std::string ref = pair.first;
            auto & cache = pair.second;
            for(size_t i = 0 ; i<cache.size(); i++ )
            {
                auto & ecctemp = cache.at(i);
                if( ecctemp.nsplit < 1 )
                {
                    if( ecctemp.ndiscordant < ndiscordant ) continue;
                }
                else
                {
                    if( ecctemp.nsplit < nsplit ||  ecctemp.score <= score ) continue;
                }
                int count_begin = cov_cache.cov_at(ref,ecctemp.begin);
                int count_end = cov_cache.cov_at(ref,ecctemp.end-1);
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
    }
    
    struct eccFinal
    {
        int begin;
        int end;
        int ndiscordant;
        int nsplit;
        float score;
        bool valid;
        float mean;
        float std;
        float begin_ratio;
        float end_ratio;
        float contigunity;
    
        void init(const EccTemp & temp)
        {
            valid = false;
            if(!temp.valid) return ;
            valid = true;
            begin = temp.begin;
            end = temp.end;
            ndiscordant = temp.ndiscordant;
            nsplit = temp.nsplit;
            score = temp.score;
            mean = -1;
            std = -1;
            begin_ratio = -1;
            end_ratio = -1;
            contigunity = -1;
        }
    
        /*
         * @abstract: fold another member of the same cluster into this one.
         *            The accumulation is the one the chained merge used when
         *            it accepted an interval.
         */
        void absorb(const EccTemp & temp)
        {
            if(!temp.valid) return ;
            if(temp.begin < begin) begin = temp.begin;
            if(temp.end > end) end = temp.end;
            if(temp.ndiscordant > ndiscordant) ndiscordant = temp.ndiscordant;
            nsplit += temp.nsplit;
            score += temp.score;
        }
        void update(const std::vector<int> & in_cov,const std::vector<int> & ext_cov)
        {
            if ( ilen <= 0 ) ilen = 100;
            if ( ext <= 0 ) ext = 200;
            if(in_cov.empty() || ext_cov.empty()) return ;
            if(in_cov.size() >1)
            {
                long long sum = 0;
                for(int x : in_cov) sum += x ;
                if( sum > 0 ) mean = double(sum)/double(in_cov.size());
                if( mean > 0 ) 
                {
                    double sd = 0 ;
                    for(int x : in_cov) sd += ((x-mean)*(x-mean));
                    std = sqrt(double(sd)/double(in_cov.size()-1));          
                }
            }
            else
            {
                mean = in_cov[0];
                std = 0;
            }

            float cov_start_in = 0;
            for(size_t i = 0 ; i<size_t(ilen) && i<in_cov.size() ; i++)
                 cov_start_in += in_cov[i];  

            float cov_end_in = 0;
            for(size_t i = in_cov.size() > size_t(ilen) ? in_cov.size()-ilen-1 : 0; i<in_cov.size() ; i++)
                 cov_end_in += in_cov[i];  

            float cov_start_ext = 0 ;
            for(size_t i = 0 ; i<size_t(ilen+ext) && i<ext_cov.size() ; i++)
                 cov_start_ext += ext_cov[i];  

            float cov_end_ext = 0 ;
            for(size_t i = ext_cov.size() > size_t(ilen+ext) ? ext_cov.size()-ilen-ext-1 : 0; i<ext_cov.size() ; i++)
                 cov_end_ext += ext_cov[i];

            if( cov_start_ext > 0.0 ) begin_ratio = cov_start_in/cov_start_ext; else begin_ratio=1.0f;
            if( cov_end_ext >0.0 ) end_ratio = cov_end_in/cov_end_ext; else end_ratio=1.0f;
            int nonzero = 0;
            for(int x : in_cov) if( x>0 ) nonzero ++;
            contigunity = 1.0 - float(nonzero)/float(in_cov.size());
        }
    };

    std::map<std::string , std::vector<eccFinal>>  eccs;
    /*
     * @abstract: group the candidates of each reference into clusters of
     *            reciprocally overlapping intervals and emit one interval per
     *            cluster.
     *
     * @notice: this used to walk the sorted candidates as a single chain,
     *          restarting the chain on every candidate that could not be added
     *          to it. That made the result depend on what sorts between two
     *          mergeable candidates: an unrelated interval in the middle reset
     *          the chain, and the candidates on either side of it were emitted
     *          separately even though they overlap each other by more than
     *          merge_factor. A spurious multi-megabase interval is the loudest
     *          way to trigger it -- on a rice eccDNA sample one such interval
     *          sitting between Chr4:1780904-1781292 (split=22) and
     *          Chr4:1780906-1781294 (split=6) kept them apart instead of
     *          merging them into Chr4:1780904-1781294 (split=28) -- but any
     *          candidate that fails the overlap test does it.
     *
     *          Clustering is independent of the order in which candidates are
     *          visited, so a candidate that merges with nothing no longer
     *          changes how its neighbours merge.
     */
    void merge()
    {
        for(auto & pair :temp_caches)
        {
            std::string ref = pair.first;
            auto & cache = pair.second;
            // invalid candidates take no part in merging, as before
            std::vector<size_t> idx;
            for(size_t i = 0 ; i<cache.size(); i++ )
                if( cache[i].valid ) idx.push_back(i);
            if( idx.empty() ) continue ;

            std::vector<size_t> parent(idx.size());
            for(size_t i = 0 ; i<parent.size(); i++ ) parent[i] = i;

            // cache is sorted by begin. Reaching 2*merge_factor reciprocal
            // overlap with a later candidate requires ov > (2f-1)*len_a and
            // ov <= len_a - (b.begin - a.begin), hence
            // b.begin - a.begin < (2-2f)*len_a. Once a later candidate starts
            // beyond that, so does every candidate after it.
            for(size_t i = 0 ; i<idx.size(); i++ )
            {
                const EccTemp & a = cache[idx[i]];
                double reach = (2.0 - 2.0*merge_factor) * double(a.end - a.begin);
                for(size_t j = i+1 ; j<idx.size(); j++ )
                {
                    const EccTemp & b = cache[idx[j]];
                    if( double(b.begin - a.begin) >= reach ) break ;
                    if( reciprocal_overlap(a,b) ) uf_unite(parent,i,j);
                }
            }

            std::map<size_t,eccFinal> clusters;
            for(size_t i = 0 ; i<idx.size(); i++ )
            {
                size_t r = uf_find(parent,i);
                auto it = clusters.find(r);
                if( it == clusters.end() )
                {
                    eccFinal one;
                    one.init(cache[idx[i]]);
                    clusters.insert(std::make_pair(r,one));
                }
                else
                    it->second.absorb(cache[idx[i]]);
            }

            std::vector<eccFinal> out;
            out.reserve(clusters.size());
            for(auto & c : clusters) out.push_back(c.second);
            std::sort(out.begin(),out.end(),
                      [](const eccFinal & x, const eccFinal & y)
                      { return x.begin < y.begin || (x.begin == y.begin && x.end < y.end); });
            for(auto & c : out) eccs[ref].push_back(c);
        }
    }

    void get_ratio()
    {
        BGIQD::MultiThread::MultiThread mt;
        for(auto & pair :eccs)
        {
            std::string ref = pair.first;
            auto & cache = pair.second;
            for(size_t i = 0 ; i<cache.size(); i++ )
            {
                
                auto & ecc = cache[i];
                eccFinal* the_ecc= &ecc;
                mt.AddJob([ref,the_ecc](){
                    int start = the_ecc->begin;
                    int end = the_ecc->end;
                    auto covs = cov_cache.cov_region(ref,start,end);
                    auto ext_covs = cov_cache.cov_region(ref,start-ext,end+ext);
                    the_ecc->update(covs,ext_covs);                 
                });
            }
        }
        mt.Start(threadnum);
        mt.WaitingStop();
    }

    void check_ratio()
    {
        for(auto & pair :eccs)
        {
            std::string ref = pair.first;
            auto & cache = pair.second;
            for(size_t i = 0 ; i<cache.size(); i++ )
            {
                auto & ecc = cache[i];
                if(ecc.begin_ratio > ratio || ecc.end_ratio > ratio )
                    ecc.valid = true;
                else
                    ecc.valid = false;
            }
        }
    }
    void print_cov(const std::string &output)
    {
        std::ofstream ofs(output);
        ofs<<"chrom\tstart\tend\tdiscordants\tsoft-clipped\tscore\tmean\tstd\tstart_ratio\tend_ratio\tcontinuity\n";
        for(auto & pair :eccs)
        {
            std::string ref = pair.first;
            auto & cache = pair.second;
            for(size_t i = 0 ; i<cache.size(); i++ )
            {
                auto & ecc = cache[i];
                if( ecc.valid ) 
                    ofs<<ref<<'\t'<<ecc.begin<<'\t'<<ecc.end<<'\t'
                       <<ecc.ndiscordant<<'\t'<<ecc.nsplit<<'\t'
                       <<ecc.score<<'\t'<<ecc.mean<<'\t'
                       <<ecc.std<<'\t'<<ecc.begin_ratio<<'\t'
                       <<ecc.end_ratio<<'\t'<<ecc.contigunity<<'\n';
            }
       }
    }
    void print_nocov(const std::string &output)
    {
        std::ofstream ofs(output);
        ofs<<"chrom\tstart\tend\tdiscordants\tsoft-clipped\tscore\n";
        for(auto & pair :eccs)
        {
            std::string ref = pair.first;
            auto & cache = pair.second;
            for(size_t i = 0 ; i<cache.size(); i++ )
            {
                auto & ecc = cache[i];
                if( ecc.valid ) 
                    ofs<<ref<<'\t'<<ecc.begin<<'\t'<<ecc.end<<'\t'
                       <<ecc.ndiscordant<<'\t'<<ecc.nsplit<<'\t'
                       <<ecc.score<<'\n';
            }
       }
    }
} ecc_cache;

void usage()
{
    std::cerr<<"\
Usage: merge_result <-c coverage.tsv>\n\
                    <-g genome.fa>\n\
                    <-i input.ecctemp.txt>\n\
                    <-o output.ecctemp.txt>\n\
                    [-N do not calculate coverage]\n\
                    [-t thread number(default 8)]\n\
                    [-n minimum number of split reads(default 0)]\n\
                    [-a allele frequency (default 0.1)]\n\
                    [-d minimum number of discordant reads(default 3)]\n\
                    [-f Merge intervals reciprocally overlapping by a fraction.(default 0.99)]\n\
                    [-e Number of bases to extend for computing the coverage ratio. Default: 200]\n\
                    [-l Number of bases inside the eccDNA breakpoint coordinates to compute the ratio. Default: 100]\n\
                    [-r Minimum in/out required coverage ratio. Default: 0.0]\n\
                    [-s minimum score (default 0.0)]"<<std::endl;
}

int main(int argc, char **argv)
{
    std::string bam;
    std::string input;
    std::string output;
    std::string genome;
    int option;
    bool NoCoverage = false;
    while((option = getopt(argc, argv, "Nc:g:c:i:o:t:n:a:d:s:hf:e:l:r:")) != -1)
    {
        switch(option)
        {
            case 'c': bam = std::string(optarg);break;
            case 'g': genome = std::string(optarg);break;
            case 'i': input = std::string(optarg);break;
            case 'o': output = std::string(optarg);break;
            case 't': threadnum = atoi(optarg);break;
            case 'n': nsplit = atoi(optarg);break;
            case 'a': af = atof(optarg);break;
            case 'f': merge_factor= atof(optarg);break;
            case 'e': ext = atoi(optarg);break;
            case 'l': ilen = atoi(optarg);break;
            case 'r': ratio = atof(optarg);break;
            case 'd': ndiscordant = atoi(optarg);break;
            case 's': score = atof(optarg);break;
            case 'N': NoCoverage = true; break;
            case 'h': usage();return 0;
            default:  usage();return 1;
        }
    }
    if( bam.empty() || input.empty() || output.empty() ||genome.empty())
    {
        usage();
        return 1;
    }
    std::cerr<<"load "<<input<<std::endl;
    ecc_cache.load(input);
    std::cerr<<"load fai of "<<genome<<std::endl;
    cov_cache.init_fai(genome);
    std::cerr<<"load "<<bam<<std::endl;
    cov_cache.init_covs(bam);
    std::cerr<<"check_1bp_cov"<<std::endl;
    ecc_cache.check_1bp_cov();   
    std::cerr<<"merge"<<std::endl;
    ecc_cache.merge();
    if( ! NoCoverage ) 
    {
        std::cerr<<"get_ratio"<<std::endl;
        ecc_cache.get_ratio();
        std::cerr<<"check_ratio"<<std::endl;
        ecc_cache.check_ratio();
    	ecc_cache.print_cov(output);
    }
    else
        ecc_cache.print_nocov(output);
    std::cerr<<"done"<<std::endl;
    return 0;
}
