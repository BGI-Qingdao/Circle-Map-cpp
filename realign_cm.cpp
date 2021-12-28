/**** std library ****/
#include <iostream>
#include <string>
#include <map>
#include <cassert>
#include <math.h>
#include <vector>
#include <set>
#include <fstream>
#include <algorithm>

/**** zlib ****/
#include <zlib.h>

/**** htslib ****/
#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "htslib/sam.h"
#include "htslib/kseq.h"
#include "htslib/regidx.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_log.h"

/**** edlib ****/
#include "utils/edlib.h"

/**** custom codes ****/
#include "utils/incr_array.h"
#include "utils/MultiThread.h"

/* ************************************************************************
 *
 * @section: utils
 *  
 *************************************************************************/

KSEQ_INIT(gzFile, gzread)

inline float round2(float var)
{
    return float( (int)(var * 100 + .5) ) / float(100) ;
}

 // prob <--> phred
static double g_qual2prob[256];
#define Phred2Prob(p) g_qual2prob[int(p)]
void InitQual2Prob()
{
    for (int i = 0; i < 256; ++i)
        g_qual2prob[i] = pow(10, -i/10.0);
}

#define bam_is_r1(b) (((b)->core.flag&BAM_FREAD1) != 0)
#define bam_is_r2(b) (((b)->core.flag&BAM_FREAD2) != 0)
#define bam_is_unmap(b) (((b)->core.flag&BAM_FUNMAP) != 0)
#define bam_is_munmap(b) (((b)->core.flag&BAM_FMUNMAP) != 0)
#define bam_is_secondary(b) (((b)->core.flag&BAM_FSECONDARY) != 0)
#define bam_is_supplementary(b) (((b)->core.flag&BAM_FSUPPLEMENTARY) != 0)


#define BLOCK_SIZE 1024*1024

/*
 * @abstract: return true if S in cigar string
 */
bool is_softclip_read(const bam1_t *b) 
{
     uint32_t* cigar_data = bam_get_cigar(b);
     for( unsigned int i = 0 ; i< (b)->core.n_cigar; i++)
         if( bam_cigar_op(cigar_data[i]) == BAM_CSOFT_CLIP )
             return true;
     return false;
}

/*
 * Get the 60 from MQ:i:60
 */
int bam_get_MQ(const bam1_t *b)
{
    uint8_t * data = NULL;
    data = bam_aux_get(b,"MQ");
    if ( data == NULL ) return 0;
    return bam_aux2i(data);
}

/*
 * @param:  bam1_t [in] alignment
 *          sc_len [out] the length of softclip
 *
 * @return: 1 for head softclip
 *          2 for tail softclip
 *          0 for no head or tail softclip
 */
inline int sc_type(const bam1_t *b, int * sc_len)
{
    if( b==NULL || sc_len == NULL) return 0;
    uint32_t* cigar_data = bam_get_cigar(b);
    if((b)->core.n_cigar == 1) return 0;
    bool begin = bam_cigar_op(cigar_data[0]) == BAM_CSOFT_CLIP;
    int  blen =  bam_cigar_oplen(cigar_data[0]);
    bool end = bam_cigar_op(cigar_data[(b)->core.n_cigar-1]) == BAM_CSOFT_CLIP ;
    int  elen =  bam_cigar_oplen(cigar_data[(b)->core.n_cigar-1]);
    if ( ! begin && ! end ) return 0;
    if ( begin && !end )
    {
        *sc_len = blen;
        return 1;
    }
    else if (!begin && end)
    {
        *sc_len = elen;
        return 2;
    }
    else
    {
        *sc_len = blen<elen?elen:blen;
        return blen<elen?2: (blen>elen?1:0);
    }
}

/*
 * @abstract: get the longest softclip only from
 *            head or tail softclip
 */
inline int longest_softcliped_base(const bam1_t *b)
{
    uint32_t* cigar_data = bam_get_cigar(b);
    assert((b)->core.n_cigar >1);
    bool begin = bam_cigar_op(cigar_data[0]) == BAM_CSOFT_CLIP;
    int  blen =  bam_cigar_oplen(cigar_data[0]);
    bool end = bam_cigar_op(cigar_data[(b)->core.n_cigar-1]) == BAM_CSOFT_CLIP ;
    int  elen =  bam_cigar_oplen(cigar_data[(b)->core.n_cigar-1]);
    assert( begin || end );
    if ( begin && !end )
        return blen;
    else if (!begin && end)
        return elen;
    else
        return blen>elen?blen:elen;
}

/*
 * @abstract: return the total length of M/X/=/I, exculde H/S/D 
 *
 * @parm: bam1_t structure and cigar in binary format.
 */
inline int infer_mapped_query_length(const bam1_t *b)
{
     int ret = 0;
     uint32_t* cigar_data = bam_get_cigar(b);
     //return bam_cigar2qlen((b)->core.n_cigar,cigar_data);
     for( unsigned int i = 0 ; i< (b)->core.n_cigar; i++)
         if( bam_cigar_op(cigar_data[i]) != BAM_CSOFT_CLIP
             && bam_cigar_op(cigar_data[i]) != BAM_CHARD_CLIP 
             && bam_cigar_op(cigar_data[i]) != BAM_CDEL)
             ret += bam_cigar_oplen(cigar_data[i]);
     return ret;
}
/*
 * @abstract: return the total length of M/X/=/D, exculde H/S/I 
 *
 * @parm: bam1_t structure and cigar in binary format.
 */
inline int infer_mapped_ref_length(const bam1_t *b)
{
     int ret = 0;
     uint32_t* cigar_data = bam_get_cigar(b);
     //return bam_cigar2rlen((b)->core.n_cigar,cigar_data);
     for( unsigned int i = 0 ; i< (b)->core.n_cigar; i++)
         if( bam_cigar_op(cigar_data[i]) != BAM_CSOFT_CLIP
             && bam_cigar_op(cigar_data[i]) != BAM_CHARD_CLIP 
             && bam_cigar_op(cigar_data[i]) != BAM_CINS)
             ret += bam_cigar_oplen(cigar_data[i]);
     return ret;
}

/*
 * @abstract: return the total length of M/X/=/D, exculde H/S/I 
 *
 * @param: cigar in char[] format like 100M4I5D35H.
 */
int infer_mapped_ref_length_from_cigar(const char * cigar)
{
    int ret = 0;
    char len_buff[4]={0,0,0,0};
    int len_index = 0 ;
    for(int i = 0;i<1000;i++)
    {
        char c = cigar[i];
        if( c>='0' && c<='9')
        {
            len_buff[len_index++]=c;
        }
        else if (c != 0)
        {
            if(c ==  'M' || c=='=' || c== 'X' || c=='D')
            {
                int op_len = atoi(len_buff);
                ret += op_len;
            }
            len_buff[0]=0;
            len_buff[1]=0;
            len_buff[2]=0;
            len_buff[3]=0;
            len_index = 0;
        }
        else if( c== 0) 
        {
           break; 
        }
    }
    return ret;
}

/*
 * @abstract: ternimate the program if bam is not sorted by reference
 */
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

inline float adaptative_myers_k(int sc_len, float edit_frac)
{
    return float(sc_len)*edit_frac;
}

inline bool non_colinearity(bool headclip, hts_pos_t pos, hts_pos_t realign_begin, hts_pos_t realign_end)
{
    if (headclip)
        return pos<realign_end /*|| pos<realign_begin*/;
    else 
        return pos>realign_begin /*|| pos>realign_end*/;
}

std::string upper(const char *line, int size)
{
    std::string ret ;
    ret.resize(size,'N');
    for(int i = 0 ; i<size ; i++)
    {
        switch(line[i])
        {
            case 'a': ret[i]='A'; break;
            case 'g': ret[i]='G'; break;
            case 'c': ret[i]='C'; break;
            case 't': ret[i]='T'; break;
            case 'A': ret[i]='A'; break;
            case 'G': ret[i]='G'; break;
            case 'C': ret[i]='C'; break;
            case 'T': ret[i]='T'; break;
            default : break;
        }
    }
    return ret;   
}

std::string seq_complementary_reverse(const std::string & line)
{
    std::string ret ;
    ret.resize(line.size(),'N');
    int index = 0;
    for( auto i = line.rbegin() ; i!= line.rend() ; i++)
    {
        if( *i == 'A' || *i == 'a' )
            ret[index++] = 'T';
        else if( *i == 'G' || *i == 'g' )
            ret[index++] = 'C';
        else if( *i == 'C' || *i == 'c' )
            ret[index++] = 'G';
        else if( *i == 'T' || *i == 't' )
            ret[index++] = 'A';
        else
            ret[index++] = 'N';
    }
    return ret;
}

void copy_hts_seq(std::string &str, const unsigned char *hts_seq, int start , int len)
{
    str.resize(len);
    for(int i = 0 ; i<len;i++)
    {
        switch(bam_seqi(hts_seq,(start+i)))
        {
            case 1: str[i]='A';break; 
            case 2: str[i]='C';break; 
            case 4: str[i]='G';break; 
            case 8: str[i]='T';break; 
            default: str[i]='N';break; 
        }
    }
}
void  copy_hts_qual(std::vector<uint8_t> & qual , const  unsigned char *qual_p,int start , int len)
{
    qual.resize(len);
    for(int i = 0 ; i<len;i++)
        qual[i] = qual_p[start+i];
}

/* ************************************************************************
 *
 * @section: global options
 *
 *************************************************************************/
//struct global_options
//{
static int threadnum = 8;
static int min_mapq = 20 ;
static int min_sc_len = 8 ;
static float edit_distance_frac= 0.05 ;
static double interval_prob_cutoff = 0.01 ;
static int nhit = 10 ;
static int gap_open = 5 ;
static int gap_extend = 1 ;
static float prob_cutoff = 0.99;
static float insert_mean = 317.87756;
static float insert_sd = 85.6517543804352;
static int std_fac = 4;
//} gopts;

/* ************************************************************************
 *
 * @section: strucures and functions to load and store discordant intervels.
 * 
 * @feature: discordant intervels are stored in ordered array
 *
 *************************************************************************/

enum MateType 
{
    LeftHook = 0,    // This read (R1) mapped at the end of eccDNA
    RightHook = 1,   // This read (R2) mapped at the begin of eccDNA
    BothHook = 2,    // Used for when merge L and R to one interval
    Unkown=999
};

struct MateInterval 
{
    MateType    type;
    hts_pos_t   reference_start;       // 0-base coordinate, include
    hts_pos_t   reference_end;         // 0-base coordinate, exclude
    hts_pos_t   next_reference_start;  // 0-base coordinate, include
    hts_pos_t   next_reference_end;    // 0-base coordinate, exclude
    double      prob; // prob = 1 - error rate
    bool        is_softclip; // true if this alignment contain soft-clip part 
                             // only non-sa mate-interval will be considered in only-discordant results
                             // incorrect r1/r2 for discordant pair also marked as true!
    // for discordant interval
    hts_pos_t   reference_end_q;       // 0-base coordinate, exclude
    hts_pos_t   next_reference_end_q;  // 0-base coordinate, exclude
                     
    std::string qname;
    // for sort algorithm
    // @NOTICE: ordered by next_reference_xxx, not reference_xxx
    bool operator <(const MateInterval & o) const 
    {
        return next_reference_start < o.next_reference_start ||
        (next_reference_start == o.next_reference_start && next_reference_end < o.next_reference_end) ; 
    }

    // for sort algorithm
    // @NOTICE: ordered by next_reference_xxx, not reference_xxx
    bool operator == (const MateInterval & o) const 
    {
        return next_reference_start == o.next_reference_start && 
               next_reference_end == o.next_reference_end ; 
    }

    // for sort algorithm
    // @NOTICE: ordered by next_reference_xxx, not reference_xxx
    bool operator > (const MateInterval & o) const 
    {
        return next_reference_start > o.next_reference_start ||
        (next_reference_start == o.next_reference_start && next_reference_end > o.next_reference_end) ; 
    }

    /*
     * create mate interval
     *
     * @param item: input alignment of one read
     *
     * Return 0 if this read is not discordant read,
     *        1 if success
     *       -1 on error 
     */
    int parse_mate_interval(const bam1_t * item)
    {
        if( !item ) return -1;
        if( bam_is_unmap(item) || bam_is_munmap(item) ) return 0;
        int mq = bam_get_MQ(item);
        //if( mq == 0 ) return 0;
        if( item->core.tid != item->core.mtid ) return 0;
        // keep only_discordant == True
        is_softclip = is_softclip_read(item);
        if( is_softclip && bam_aux_get(item,"SA") != NULL ) return 0;
        // check discordant reads (R2F1 orientation, R2 read)
        if( bam_is_rev(item) && ! bam_is_mrev(item) && item->core.pos < item->core.mpos )
        {
            type = MateType::RightHook;
            reference_start = item->core.pos;
            reference_end = reference_start + infer_mapped_ref_length(item);
            int infer_read_len = infer_mapped_query_length(item);
            reference_end_q = reference_start + infer_read_len;
            next_reference_start = item->core.mpos;
            next_reference_end = next_reference_start + item->core.l_qseq;
            next_reference_end_q = next_reference_start + infer_read_len;
            prob = 1.0 - Phred2Prob(mq);
            qname=std::string(bam_get_qname(item));
            if(bam_is_r1(item)) is_softclip = true;
            //@debug
            //std::cerr<<qname<<'\t'<<reference_start<<'\t'<<reference_end<<'\t'
            //            <<next_reference_start<<'\t'<<next_reference_end<<'\n';
            return 1;
        }
        // R2F1 when iterating trough R1 read
        if( ! bam_is_rev(item) &&  bam_is_mrev(item) && item->core.pos > item->core.mpos )
        {
            type = MateType::LeftHook;
            reference_start = item->core.pos;
            reference_end = reference_start + infer_mapped_ref_length(item);
            int infer_read_len = infer_mapped_query_length(item);
            reference_end_q = reference_start + infer_read_len;
            next_reference_start = item->core.mpos;
            next_reference_end = next_reference_start +  item->core.l_qseq;
            next_reference_end_q = next_reference_start + infer_read_len;
            prob = 1.0 - Phred2Prob(mq);
            qname=std::string(bam_get_qname(item));
            if(bam_is_r2(item)) is_softclip = true;
            //@debug
            //std::cerr<<qname<<'\t'<<reference_start<<'\t'<<reference_end<<'\t'
            //            <<next_reference_start<<'\t'<<next_reference_end<<'\n';
            return 1;
        }
        return 0;
    }
};


/* ************************************************************************
 *
 * @section: strucures and functions to load and store eccDNA datas.
 *
 *************************************************************************/

struct eccTemp
{
    std::string qname;
    hts_pos_t   pos_begin;
    hts_pos_t   pos_end;  
    hts_pos_t   leftmost_sa;
    hts_pos_t   cut_begin; // 0base include 
    hts_pos_t   cut_end;   // 0base exclude
    float       score;// score = round(xxx,2)

    // for sort algorithm
    // @NOTICE: ordered by cut_xxx 
    bool operator <(const eccTemp & o) const 
    {
        return cut_begin < o.cut_begin ||
        (cut_begin == o.cut_begin && cut_end < o.cut_end) ; 
    }

    // for sort algorithm
    // @NOTICE: ordered by cut_xxx
    bool operator == (const eccTemp & o) const 
    {
        return cut_begin == o.cut_begin && 
               cut_end == o.cut_end; 
    }

    // for sort algorithm
    // @NOTICE: ordered by cut_xxx
    bool operator > (const eccTemp& o) const 
    {
        return cut_begin > o.cut_begin ||
        ( cut_begin == o.cut_end &&  cut_end> o.cut_end ) ; 
    }
    /*
     * create ecc proof
     *
     * @param item: input alignment of one read
     *
     * Return 0 if this read is not split reads,
     *        1 if success
     *       -1 on error 
     */
    int parse_ecc( const bam1_t * item, bam_hdr_t* header )
    {
        if( !item ) return -1;
        uint8_t * data =  bam_aux_get(item,"SA");
        if( !is_softclip_read(item) || data == NULL ) return 0;
        kstring_t sa_tmp = { 0, 0, NULL };
        kputs( bam_aux2Z(data),&sa_tmp );
        int *offset = NULL;
        int num_item = 0;
        // SA:Z:chr8,66989,-,89M61S,0,0;
        if ((offset = ksplit(&sa_tmp, ',', &num_item)) == NULL) return -1;
        if (num_item <5) return -1;
        char sa_orientation = *(sa_tmp.s + offset[2]);
        if( sa_orientation == '+' &&  bam_is_rev(item) ) return 0;
        if( sa_orientation == '-' && !bam_is_rev(item) ) return 0;
        if( item->core.tid != sam_hdr_name2tid(header,sa_tmp.s + offset[0]) )
            return 0;
        int sq = atoi(sa_tmp.s + offset[4]);
        //@NOTICE: in python code, here is sa > min_mapq
        if ( sq <= min_mapq )
        {
            free(sa_tmp.s);
            free(offset);
            return 0;
        }
        qname=std::string(bam_get_qname(item));
        pos_begin =  item->core.pos;
        pos_end = pos_begin + infer_mapped_ref_length(item);
        leftmost_sa = atol(sa_tmp.s + offset[1]); // 1base
        score = round2( float( longest_softcliped_base(item) ) *(1-Phred2Prob(sq)) );

        if( pos_begin < leftmost_sa )
        {
            cut_begin = pos_begin;
            const char * cigar = sa_tmp.s + offset[3];
            cut_end = leftmost_sa -1 + infer_mapped_ref_length_from_cigar(cigar) -1 ;
        }
        else if( pos_begin > leftmost_sa )
        {
            cut_begin = leftmost_sa -1 ; // 0base
            cut_end = pos_end; 
        }
        free(sa_tmp.s);
        free(offset);
        return 1;
    }
};


/* ************************************************************************
 *
 * @section: strucures and functions to handle Soft-Clip reads.
 *
 *************************************************************************/

struct scRead
{
    struct sc_part
    {
        bool headclip ;
        std::string seqs;
        std::vector<uint8_t> qual;
    } sc_info;
    std::string qname;
    bool is_reverse;
    hts_pos_t pos_begin;
    hts_pos_t pos_end;
    /*
     * create Softclip reads
     *
     * @param item: input alignment of one read
     *
     * Return 0 if this read is not softclip reads,
     *        1 if success
     *       -1 on error 
     */
    int parse_sc(const bam1_t * item, const bam_hdr_t* header,int min_sc_len)
    {
        if( !item ) return -1;
        int sc_len = 0;
        int sctp = sc_type(item,&sc_len);
        uint8_t * data =  bam_aux_get(item,"SA");
        if( sctp == 0 || data != NULL ) return 0;
        if(sc_len < min_sc_len ) return 0;
        // A true softclip reads here
        pos_begin = item->core.pos;
        pos_end = pos_begin + infer_mapped_ref_length(item);
        is_reverse = bam_is_rev(item);
        qname=std::string(bam_get_qname(item));
        assert(sc_len>0);
        if( sctp == 1 )
        {
            sc_info.headclip = true;
            copy_hts_seq(sc_info.seqs ,bam_get_seq(item),0,sc_len);
            copy_hts_qual(sc_info.qual,bam_get_qual(item),0,sc_len);
        } 
        else
        {
            sc_info.headclip = false;
            copy_hts_seq(sc_info.seqs ,bam_get_seq(item),(item)->core.l_qseq-sc_len,sc_len);
            copy_hts_qual(sc_info.qual,bam_get_qual(item),(item)->core.l_qseq-sc_len,sc_len);
        } 
        return 1;
    }
};

/* ************************************************************************
 *
 * @section: global caching structure.
 *
 *************************************************************************/

typedef BGIQD::INCRARRAY::IncrArray<MateInterval> MateArray;
typedef BGIQD::INCRARRAY::IncrArray<eccTemp>      eccArray;
typedef BGIQD::INCRARRAY::IncrArray<scRead>       scArray;
struct all_caches
{
    std::map<std::string, MateArray> chromesomes_mates;
    std::map<std::string, eccArray>  chromesomes_eccs;
    std::map<std::string, scArray>   chromesomes_scs;
    std::map<std::string, kstring_t>    chromesomes;   
    private:
    // add mate interval into cache
    inline void add_mate(const std::string & chrname, const MateInterval & tmp) 
    {
        auto itr = chromesomes_mates.find(chrname);
        if( itr == chromesomes_mates.end() )
        {
            chromesomes_mates[chrname].Init(BLOCK_SIZE);
            chromesomes_mates[chrname].push_back(tmp);
        }
        else
            itr->second.push_back(tmp);
    }
    inline void add_ecc(const std::string & chrname, const eccTemp& tmp) 
    {
        auto itr = chromesomes_eccs.find(chrname);
        if( itr == chromesomes_eccs.end() )
        {
            chromesomes_eccs[chrname].Init(BLOCK_SIZE);
            chromesomes_eccs[chrname].push_back(tmp);
        }
        else
            itr->second.push_back(tmp);
    }
    inline void add_sc(const std::string & chrname, const scRead & tmp) 
    {
        auto itr = chromesomes_scs.find(chrname);
        if( itr == chromesomes_scs.end() )
        {
            chromesomes_scs[chrname].Init(BLOCK_SIZE);
            chromesomes_scs[chrname].push_back(tmp);
        }
        else
            itr->second.push_back(tmp);
    }

    inline void add_chr(const kseq_t *ks)
    {
        chromesomes[std::string(ks->name.s)] = {0,0,NULL};
        kputsn(ks->seq.s,ks->seq.l, &(chromesomes[std::string(ks->name.s)]));
    }
        
    public:

    void free_mates()
    {
        for(auto & chr : chromesomes_mates)
            chr.second.deep_clean();
    }

    void free_eccs()
    {
        for(auto & chr : chromesomes_eccs)
            chr.second.deep_clean();
    }

    void free_scs()
    {
        for(auto & chr : chromesomes_scs)
            chr.second.deep_clean();
    }

    void load_genome(const char * filename)
    {
        gzFile fp_fa = gzopen(filename, "r");
        kseq_t *ks;
        ks = kseq_init(fp_fa);
        while ( kseq_read(ks) >= 0) 
            add_chr(ks);
        kseq_destroy(ks);
        gzclose(fp_fa);
    }
    /*
     * @abstract: get sequence slice from interval chrname -- [start,end)
     *
     * @notice: all bases are capitalized
     *
     */ 
    std::string get_genome_slice(const std::string & chrname, hts_pos_t start, hts_pos_t end)
    {
         const auto & itr =  chromesomes.find(chrname);
         if( itr == chromesomes.end() )
         {
              hts_log_error("Error: genome not match bam!  \"%s\" not found. exit... ", chrname.c_str());
              return "";
         }
         if (end<= start) return "";
         if(start < 0 ) return "";
         if(end >= hts_pos_t(itr->second.l)) end = hts_pos_t(itr->second.l)-1; 
         return upper(itr->second.s + start , end-start);
    }
    
    //
    // @abstract: load candidate eccDNA.bam and cache 
    //                 a) mate intervals from discordant reads pair
    //                 b) soft-clip reads from SA reads.
    //                 c) ecc intervals from split reads
    //
    // @notice: input bam must be sorted by mapping position.
    //
    // @param: sort_candidates_bam_filename -- the input bam filename
    //
    void load_mates_sa_sr(const std::string & sort_candidates_bam_filename)
    {
        htsThreadPool p = {NULL,0};
        // try to open bam file
        samFile *in = 0;
        if ((in = sam_open(sort_candidates_bam_filename.c_str(), "r")) == 0) 
        {
             hts_log_error("Failed to open sorted bam : \"%s\" for reading", sort_candidates_bam_filename.c_str());
             exit(1);
        }
        // prepare thread pool
        if( threadnum > 1 )
        {
            if (!(p.pool = hts_tpool_init(threadnum ) ) ) 
            {
                hts_log_error("Create tpool failed");
                exit(1);
            }                                                   
            hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        }
        // handle alignment one by one
        bam1_t *b = bam_init1();
        bam_hdr_t * header = NULL;
        header = sam_hdr_read(in);
        if (! header ) 
        {
            hts_log_error("Failed to load header from : \"%s\"", sort_candidates_bam_filename.c_str());
            exit(1);
        }
        guarentee_bam_format(header);
        int r;
        MateInterval temp;
        eccTemp etemp;
        scRead stemp;
        while ((r = sam_read1(in, header, b)) >= 0) 
        { // read one alignment from `in'
            if( b->core.qual >= min_mapq )
            {
                if (temp.parse_mate_interval(b) >0 ) // handle discordant reads
                    add_mate(std::string(sam_hdr_tid2name(header,b->core.tid)),temp);
                if ( etemp.parse_ecc(b,header)>0 )
                    add_ecc(std::string(sam_hdr_tid2name(header,b->core.tid)),etemp);
                if ( stemp.parse_sc(b,header,min_sc_len)>0 )
                    add_sc(std::string(sam_hdr_tid2name(header,b->core.tid)),stemp);
            }
        }
        // release resources
        sam_hdr_destroy(header);
        bam_destroy1(b);
        sam_close(in);
        if( threadnum > 1 )
            hts_tpool_destroy(p.pool);
        // handle sam_read1 return after release resources
        if (r < -1) 
        {
            hts_log_error("Error reading file : \"%s\"", sort_candidates_bam_filename.c_str());
            exit(1);
        }
        // print log
        long long total = 0;
        for( const auto &itr : chromesomes_mates)
            total += itr.second.size();
        hts_log_info("total mate intervals : %llu ",total);
        return ;
    }
    
} caches;

/* ************************************************************************
 *
 * @section: strucures and functions to store and manage peaks.
 *
 * @feature: peaks are stored in a ordered array.
 *
 *************************************************************************/

struct eccFinal
{
     hts_pos_t cut_begin;
     hts_pos_t cut_end;
     int       split_read_count;
     int       discordant_read_count;
     float     score;
     std::set<std::string> s_qnames;
     void init(const eccTemp & temp)
     {
         cut_begin = temp.cut_begin;
         cut_end = temp.cut_end;
         score = temp.score;
         //assert(temp.score>0);
         discordant_read_count = 0;
         split_read_count = 1;
         s_qnames.clear();
         s_qnames.insert(temp.qname);
     }
     /*
      * @abstract: check if temp is exact same with this
      * 
      * @notice: please add ordered data
      * 
      * @return: true if add succ due to temp equal with this.
      *          false otherwise.
      */  
     bool try_add_next(const eccTemp & temp)
     {
         //check overlap
         if ( temp.cut_begin != cut_begin
           || temp.cut_end != cut_end ) return false;
         split_read_count ++ ;
         s_qnames.insert(temp.qname);
         score += temp.score;
         return true;
     }
};

struct discordant
{
     hts_pos_t cut_begin;
     hts_pos_t cut_end;
     int       discordant_read_count;
     std::set<std::string> qnames;
     void init(const MateInterval & temp)
     {
          discordant_read_count = 0;
          qnames.clear();
          if( temp.is_softclip ) return;
          discordant_read_count = 1;
          qnames.insert(temp.qname);
          if( temp.reference_start < temp.next_reference_start )
          {
              cut_begin = temp.reference_start;
              cut_end = temp.next_reference_end_q;
          }
          else if (temp.reference_start >temp.next_reference_start)
          {
              cut_begin = temp.next_reference_start;
              cut_end = temp.reference_end_q;
          }
          else
          {
              discordant_read_count=0;
              qnames.clear();
          }
     }

     /*
      * @abstract: check if temp is overlap with this
      * 
      * @notice: please add ordered data
      * 
      * @return: true if add succ due to temp overlap with this.
      *          false otherwise.
      */  
     bool try_add_next(const MateInterval & temp)
     {
          if( temp.is_softclip ) return true ; // skip all softclip pair
          if( discordant_read_count == 0 ) 
          {
              init(temp);
              return true;
          }
          if( temp.reference_start < temp.next_reference_start ) 
          {
              if( temp.reference_start < cut_end )
              {
                  discordant_read_count ++ ;
                  if ( temp.next_reference_end_q > cut_end )
                      cut_end = temp.next_reference_end_q;
                  qnames.insert(temp.qname);
                  return true;
              }
          }
          else if ( temp.reference_start > temp.next_reference_start )
          {
              if( temp.next_reference_start < cut_end )
              {
                  discordant_read_count ++ ;
                  if ( temp.reference_end_q > cut_end )
                      cut_end = temp.reference_end_q;
                  qnames.insert(temp.qname);
                  return true;
              }
          }
          return false;
     }
};

/*
 * @abstract:  details of one candidate interval.
 *             core structure of this realign algorithm.
 *
 */
struct ReAlignItem
{
    public:
    ReAlignItem(){}
    ReAlignItem(hts_pos_t b,hts_pos_t e):begin(b),end(e){}
    hts_pos_t begin;    //0 base. The same as `bedtools genomecov -bg`
    hts_pos_t end;      //0 base. The same as `bedtools genomecov -bg`
    /***********************************************************************************/
    std::vector<MateInterval>  mate_intervals;   // mate interval by discordant reads
    std::vector<MateInterval>  realign_interval; // realign interval that merged and filted by mate intevals
    std::vector<scRead>        scs;              // softclip reads that waiting to realign
    std::vector<eccFinal>      final_eccs;       // final ecc intervals
    std::vector<discordant>    only_discordant;  // discordant interval if no final eccs found.
    std::vector<eccTemp>       eccs;             // intermidate ecc intervals
    /***********************************************************************************/
    int realign_interval_valid;
    private:
    MateType all_type ;
    
    public:
    // for sort algorithm
    bool operator <(const ReAlignItem & o) const { return begin < o.begin || (begin==o.begin && end < o.end) ; }
    bool operator == (const ReAlignItem & o) const { return begin == o.begin && end == o.end ; }
    bool operator > (const ReAlignItem & o) const { return begin > o.begin || (begin==o.begin && end > o.end) ; }

     /*
      * @abstract: set realign interval by mate_intervals
      *
      * @pram: interval_prob_cutoff -- only consider realign interval if prob > interval_prob_cutoff 
      * 
      * @return: 0 if no >99% probility interval,
      *          1 otherwise
      *
      * @warn: this function should NOT be called twice.
      */
    int set_realign_interval(double interval_prob_cutoff) 
    {
        realign_interval_valid = 0;
        if ( mate_intervals.empty() )
            return 0;
        if ( mate_intervals.size() == 1 )
        {
            if( mate_intervals[0].prob == 0 ) return 0;
            realign_interval_valid = 1;
            realign_interval.push_back(mate_intervals[0]);
            realign_interval[0].prob = 1.0;
            all_type = realign_interval[0].type;
            
            return 1;
        }
        // sort by mate_reference_pos
        std::sort(mate_intervals.begin(), mate_intervals.end());
        // merge overlaped intervals and sum their probabilites
        std::vector<MateInterval> candidate_mate_intervals;
        MateInterval temp = mate_intervals.at(0);
        all_type = temp.type;
	for(size_t i = 1; i < mate_intervals.size(); i++)
        {
            const auto & curr = mate_intervals.at(i);
            if ( curr.type != all_type ) all_type =  MateType::BothHook;
            //if ( curr.next_reference_start >= temp.next_reference_end ) // new candidate mate interval
            if ( curr.next_reference_start > temp.next_reference_end ) // new candidate mate interval
            {
                candidate_mate_intervals.push_back(temp);
                temp = curr ;
            }
            else 
            { 
                // merge into previous mate interval
                temp.prob += curr.prob;
                //the python code use last instead of max
                if( curr.next_reference_end > temp.next_reference_end )
                    temp.next_reference_end = curr.next_reference_end;
                if( curr.type != temp.type )
                    temp.type = MateType::BothHook;
            }
        }
        candidate_mate_intervals.push_back(temp);
        double total_prob = 0.0f;
        for( const auto & mi : candidate_mate_intervals)
            total_prob += mi.prob;
        if( total_prob == 0 ) return 0;
        for( auto & mi : candidate_mate_intervals)
            mi.prob = float(mi.prob)/float(total_prob);
        int ret = 0;
        for( const auto & mi : candidate_mate_intervals)
        {
            if ( mi.prob >= double(interval_prob_cutoff))
            {
                realign_interval_valid = 1;
                realign_interval.push_back(mi);
                ret = 1;
            }
        }   
        return ret;
    }
    //
    // @abstract: expand the realign interval
    //
    // @param: interval_extand_size -- Eg: Mean(insert_size) * 3 SD(insert_size)
    //
    void expand_realign_interval(int interval_extand_size) 
    {
        assert(realign_interval_valid == 1 );
        for( size_t i = 0; i<realign_interval.size() ; i++ )
        {
            hts_pos_t start = realign_interval[i].next_reference_start > interval_extand_size ? 
                              realign_interval[i].next_reference_start-interval_extand_size :0 ;
            hts_pos_t end = realign_interval[i].next_reference_end+interval_extand_size;
            if ( all_type == MateType::LeftHook || all_type == MateType::BothHook)
            	realign_interval[i].next_reference_start = start;
            if ( all_type == MateType::RightHook || all_type == MateType::BothHook)
                realign_interval[i].next_reference_end = end;
            //if (realign_interval[i].type == MateType::LeftHook || realign_interval[i].type == MateType::BothHook) 
            //    realign_interval[i].next_reference_start = start;
            //if (realign_interval[i].type == MateType::RightHook ||realign_interval[i].type == MateType::BothHook)
            //    realign_interval[i].next_reference_end = end;
        }
        // remerge interval
        if( realign_interval.size() > 1 ) 
        {
            std::vector<MateInterval> temp;
            MateInterval prev = realign_interval[0];
            for( size_t i = 1; i<realign_interval.size() ; i++ )
            {
                const auto & next = realign_interval.at(i);
                if( next.next_reference_start  > prev.next_reference_end )
                {
                    temp.push_back(prev);
                    prev = next ;
                }
                else
                {
                    if(prev.next_reference_end <next.next_reference_end)
                        prev.next_reference_end = next.next_reference_end;
                }
            }
            temp.push_back(prev);
            realign_interval.swap(temp);
        }
    }

    // add mate interval into mate_intervals
    void hook_mate(const MateInterval &c) 
    {
        mate_intervals.push_back(c);
        return;
    }

    // add eccs if valid
    void hook_ecc(const eccTemp & e, bool force)
    {
        ////@debug
        //8260985, 8262149
        //6512791 6513277
        ////@debug
        //if( temp.cut_begin == 8664950 && temp.cut_end == 8665616)
        //    std::cerr<<temp.qname<<'\t'<<temp.score<<std::endl;
        if( force )
        {
            eccs.push_back(e);
            return;
        }
        // ignore this is no valid realign(mate) interval
        if( realign_interval_valid == 0 ) return ;
        for( size_t i = 0; i<realign_interval.size() ; i++ )
        {
            const auto & realign = realign_interval.at(i);
            // check if fit realign  inteval
            if( realign.next_reference_start < e.leftmost_sa && e.leftmost_sa <  realign.next_reference_end )
            {
                //if( e.cut_begin == 6512791&& e.cut_end == 6513277 )
                //{
                //    std::cerr<<e.qname<<std::endl;
                //    std::cerr<<realign.next_reference_start<<std::endl;
                //    std::cerr<<realign.next_reference_end<<std::endl;
                //    std::cerr<<begin<<std::endl;
                //    std::cerr<<end<<std::endl;
                //}
                eccs.push_back(e);
                //return;
            }
        }
    }
    // add sc if valid
    void hook_sc(const scRead & s)
    {
        // ignore this is no valid realign(mate) interval
        if( realign_interval_valid == 0 ) return ;
        for( size_t i = 0; i<realign_interval.size() ; i++ )
        {
            const auto & realign = realign_interval.at(i);
            // check if fit realign  inteval
            if( non_colinearity(s.sc_info.headclip, s.pos_begin, 
                                realign.next_reference_start,
                                realign.next_reference_end ) )
            {
                scs.push_back(s);
                return;
            }
        }
    }
    bool has_eccs() const
    {
        if( realign_interval_valid == 0 ) return false;
        return true;
    }

    void assign_discordant()
    {
        float max_dist = float(insert_mean) + float(5.0 * insert_sd);
        if( eccs.empty() ) return ;
        //sort eccs
        if( eccs.size() > 1 )
            std::sort(eccs.begin(),eccs.end());
        //merge eccs to final_eccs
        eccFinal temp;
        temp.init(eccs[0]);
        for( size_t i = 1; i<eccs.size();i++)
        {
            if(false == temp.try_add_next(eccs[i]))
            {
                final_eccs.push_back(temp);
                temp.init(eccs[i]);
            }
        }
        final_eccs.push_back(temp);
        //assign discordant reads
        for(auto & fecc : final_eccs)
        {
            //@debug
            ///if( fecc.cut_begin == 8850850 && fecc.cut_end == 8851833)
            ///{
            ///    for( const auto & x : fecc.s_qnames)
            ///        std::cerr<<x<<std::endl;
            ///}
            std::set<std::string> d_qnames;
            for(size_t i = 0; i<mate_intervals.size(); i++)
            {
                const auto &mi = mate_intervals.at(i);
                //if( mi.qname == "V350014168L4C002R0470929844")
                //     std::cout<<1<<std::endl;
                if( mi.is_softclip ) continue;
                //hts
                if( mi.reference_start < mi.next_reference_start ) 
                {
                    if( mi.reference_start >fecc.cut_begin
                     && mi.reference_start-fecc.cut_begin<max_dist
                     && mi.next_reference_end_q < fecc.cut_end
                     && fecc.cut_end - mi.next_reference_end_q<max_dist)
                    {
                        //fecc.discordant_read_count++;
                        d_qnames.insert(mi.qname);
                        // @debug
                        //if( fecc.cut_begin == 7918273 && fecc.cut_end == 7918655 )
                        //    std::cerr<<mi.qname<<std::endl;
                    }
                   
                }
                else
                {
                    if( mi.next_reference_start >fecc.cut_begin
                     && mi.next_reference_start-fecc.cut_begin<max_dist
                     && mi.reference_end_q < fecc.cut_end
                     && fecc.cut_end - mi.reference_end_q<max_dist)
                    {
                        //fecc.discordant_read_count++;
                        d_qnames.insert(mi.qname);
                        // @debug
                        //if( fecc.cut_begin == 7918273 && fecc.cut_end == 7918655 )
                        //    std::cerr<<mi.qname<<std::endl;
                    }
                }
            }
            fecc.discordant_read_count = d_qnames.size();
            //debug
            //if ( fecc.cut_begin == 15104605 && fecc.cut_end == 15106368 )
            //    for(const auto x :d_qnames)std::cerr<<x<<std::endl;
        }
    }
    
    void discordant_only()
    {
        discordant temp;
        temp.init(mate_intervals[0]); 
        for( size_t i = 1; i<mate_intervals.size(); i++)
        {
            if( false == temp.try_add_next(mate_intervals[i]) )
                only_discordant.push_back(temp);
        }
        if(temp.discordant_read_count >0)
            only_discordant.push_back(temp);
    }
};

typedef ReAlignItem region_position;

typedef std::vector<region_position>  region_of_chromesome;

struct all_region
{
    std::map<std::string, region_of_chromesome> chromesomes;
    
    //
    // @abstract: load chunked candidate intervals from peak.bed and sort them by coordinate
    // 
    // @param: peak_bed_filename -- the filename of peak.bed
    // 
    void load_peaks(const std::string & peak_bed_filename)
    {
        // load peaks.bed
        regidx_t *idx = regidx_init(peak_bed_filename.c_str(),NULL,NULL,0,NULL);
        regitr_t *itr = regitr_init(idx);
        // create chunks of candidate intervals
        unsigned long long total = 0;
        while ( regitr_loop(itr) )
        {
            total++;
            std::string chrname = std::string(itr->seq);
            hts_pos_t beg = itr->beg;
            hts_pos_t end = itr->end+1; // regidx use 0 base in end. Keep the same as peak.bed.
            if ( end - beg <= 500 )
                chromesomes[chrname].push_back(ReAlignItem(beg,end)); 
            else 
            {
                while( beg < end )
                {
                    hts_pos_t end300 = beg + 300 ;
                    chromesomes[chrname].push_back(ReAlignItem(beg,end300));
                    beg += 300 ;
                }
            }
        }
        hts_log_info("peaks.bed intervals : %llu" ,total);
        regidx_destroy(idx);                                                    
        regitr_destroy(itr);                                                    
        total = 0;
        for( auto & itr : chromesomes) 
        { 
            std::sort(itr.second.begin(),itr.second.end());
            total += itr.second.size();
        }
        hts_log_info("chunked candidate intervals : %llu ",total);
        return ;
    }

    //
    // @abstract: merge and filter mate intervals to get realign interval (>99% prob)
    // 
    void get_realign_intervals()
    {
        for( auto & itr : chromesomes)
        {
             region_of_chromesome & chr = itr.second;
             for( size_t i = 0 ; i< chr.size(); i++)
             {
                 if( chr[i].set_realign_interval(interval_prob_cutoff) == 1)
                 {
                     chr[i].expand_realign_interval(insert_mean+insert_sd*std_fac);
                 }
             }
        }
    }
    //
    // @abstract: print all valid candidate intervals and their cooresponded 
    //            realign interval for debuging
    //
    void debug_print_realign_intervals()
    {
        for( const auto & itr : chromesomes)
        {
             const auto & chrname = itr.first;
             const region_of_chromesome & chr = itr.second;
             for ( size_t i = 0 ; i< chr.size(); i++ )
             {
                 const ReAlignItem & item = chr.at(i);
                 if( item.realign_interval_valid == 1)
                 {
                     for( size_t j = 0 ; j < item.realign_interval.size(); j ++)
                         std::cerr<<chrname<<'\t'<<item.begin<<'\t'<<item.end<<'\t'
                              <<item.realign_interval.at(j).next_reference_start<<'\t'
                              <<item.realign_interval.at(j).next_reference_end<<'\n';
                 }
             }
        }
    }
} peaks;

/* ************************************************************************
 *
 * @section:  functions to merge caches to peaks.
 *
 *************************************************************************/

void merge_candidate_mate()
{
    for(const auto & itr : caches.chromesomes_mates) 
    {
        const std::string & chrname = itr.first;
        const MateArray & mate_intervals = itr.second;
        if( peaks.chromesomes.find(chrname) == peaks.chromesomes.end())
            continue;
        region_of_chromesome & candidates =  peaks.chromesomes.at(chrname);
        size_t itr_m = 0, itr_c = 0;
        for( ; itr_m < mate_intervals.size() ; itr_m++)
        {
            const auto & mate = mate_intervals.at(itr_m);
            hts_pos_t m_begin = mate.reference_start, m_end = mate.reference_end;
            hts_pos_t c_begin = 0, c_end = 0;
            bool mate_is_left = false;
            while(1)
            {
                if(itr_c >= candidates.size())
                    break;
                const auto & candidate = candidates.at(itr_c);
                c_begin = candidate.begin;
                c_end = candidate.end;
                //std::tie(c_begin,c_end) = candidates.at(itr_c);
                if( m_end <= c_begin ) 
                {
                    mate_is_left = true;
                    break;
                } 
                if ( m_begin < c_end ) // overlap now
                    break ;
                itr_c ++ ;
            }
            if(itr_c >= candidates.size())
                break;
            if( mate_is_left )
                continue;
            // overlap now
            auto & candidate = candidates.at(itr_c);
            //@debug
            //std::cerr<<mate.qname<<'\t'<<mate.reference_start<<'\t'
            //                           <<mate.reference_end<<'\t'
            //                           <<mate.next_reference_start<<'\t'
            //                           <<mate.next_reference_end<<'\t'
            //                           <<std::endl;
            candidate.hook_mate(mate);
            // check if this mate interval also overlap with next candidate
            if ( itr_c < candidates.size()-1 ) 
            {
                auto & next_candidate  = candidates.at(itr_c+1);
                c_begin = next_candidate.begin;
                c_end = next_candidate.end;
                if( m_end > c_begin && m_begin < c_end )
                {
                    next_candidate.hook_mate(mate);
                }
            }
        }
    }
    caches.free_mates();
}

void merge_candidate_split_reads()
{
    for(const auto & itr : caches.chromesomes_eccs) 
    {
        const std::string & chrname = itr.first;
        const eccArray & eccs = itr.second;
        if( peaks.chromesomes.find(chrname) == peaks.chromesomes.end())
            continue;
        region_of_chromesome & candidates =  peaks.chromesomes.at(chrname);
        size_t itr_e = 0, itr_c = 0;
        for( ; itr_e < eccs.size() ; itr_e++)
        {
            const auto & ecc = eccs.at(itr_e);
            hts_pos_t e_begin = ecc.pos_begin , e_end = ecc.pos_end;
            hts_pos_t c_begin = 0, c_end = 0;
            bool ecc_is_left = false;
            while(1)
            {
                if(itr_c >= candidates.size())
                    break;
                const auto & candidate = candidates.at(itr_c);
                c_begin = candidate.begin;
                c_end = candidate.end;
                //std::tie(c_begin,c_end) = candidates.at(itr_c);
                if( e_end <= c_begin ) 
                {
                    ecc_is_left = true;
                    break;
                } 
                if ( e_begin < c_end ) // overlap now
                    break ;
                itr_c ++ ;
            }
            if(itr_c >= candidates.size())
                break;
            if( ecc_is_left )
                continue;
            auto & candidate = candidates.at(itr_c);
            candidate.hook_ecc(ecc,false);
            if ( itr_c < candidates.size()-1 ) 
            {
                auto & next_candidate  = candidates.at(itr_c+1);
                c_begin = next_candidate.begin;
                c_end = next_candidate.end;
                if( e_end > c_begin && e_begin < c_end )
                {
                    next_candidate.hook_ecc(ecc,false);
                }
            }
        }
    }
    caches.free_eccs();
}

void merge_candidate_sc()
{
    for(const auto & itr : caches.chromesomes_scs)
    {
        const std::string & chrname = itr.first;
        const scArray & scs = itr.second;
        if( peaks.chromesomes.find(chrname) == peaks.chromesomes.end())
            continue;
        region_of_chromesome & candidates =  peaks.chromesomes.at(chrname);
        size_t itr_s = 0, itr_c = 0;
        for( ; itr_s < scs.size() ; itr_s++)
        {
            const auto & sc = scs.at(itr_s);
            hts_pos_t s_begin = sc.pos_begin , s_end = sc.pos_end;
            hts_pos_t c_begin = 0, c_end = 0;
            bool sc_is_left = false;
            while(1)
            {
                if(itr_c >= candidates.size())
                    break;
                const auto & candidate = candidates.at(itr_c);
                c_begin = candidate.begin;
                c_end = candidate.end;
                //std::tie(c_begin,c_end) = candidates.at(itr_c);
                if( s_end <= c_begin ) 
                {
                    sc_is_left = true;
                    break;
                } 
                if ( s_begin < c_end ) // overlap now
                    break ;
                itr_c ++ ;
            }
            if(itr_c >= candidates.size())
                break;
            if( sc_is_left )
                continue;
            auto & candidate = candidates.at(itr_c);
            candidate.hook_sc(sc);
            if ( itr_c < candidates.size()-1 ) 
            {
                auto & next_candidate  = candidates.at(itr_c+1);
                c_begin = next_candidate.begin;
                c_end = next_candidate.end;
                if( s_end > c_begin && s_begin < c_end )
                {
                    next_candidate.hook_sc(sc);
                }
            }
        }
    }
    caches.free_scs();
    return ;
}


/* ************************************************************************
 *
 * @section:  pssm module
 *
 *************************************************************************/

struct AGCTFreq 
{
    float freq[4]; // AGCT freq
    void init(const std::string & seq)
    {
        float total = 0, A=0,G=0,C=0,T=0;
        for(char c : seq)
        {
            switch(c)
            {
                case 'A': A++; total++;break;
                case 'G': G++; total++;break;
                case 'C': C++; total++;break;
                case 'T': T++; total++;break;
                default:break;
            }
        }
        if(total == 0)
        {
            freq[0]=0;freq[1]=0;freq[2]=0;freq[3]=0;
        }
        else
        {
            freq[0]=A/total;freq[1]=G/total;freq[2]=C/total;freq[3]=T/total;
        }
    }
};

inline float base_prob_match(char base,char qual,const AGCTFreq & used_ref_freq)
{
    switch(base)
    {
        case 'A': return log2((1-Phred2Prob(qual))/used_ref_freq.freq[0]);
        case 'G': return log2((1-Phred2Prob(qual))/used_ref_freq.freq[1]);
        case 'C': return log2((1-Phred2Prob(qual))/used_ref_freq.freq[2]);
        case 'T': return log2((1-Phred2Prob(qual))/used_ref_freq.freq[3]);
        default : return 0;
    }
    return 0;
}

inline float base_prob_mismatch(char base,char qual,const AGCTFreq & used_ref_freq)
{
    switch(base)
    {
        case 'A': return log2((Phred2Prob(qual)/3)/used_ref_freq.freq[0]);
        case 'G': return log2((Phred2Prob(qual)/3)/used_ref_freq.freq[1]);
        case 'C': return log2((Phred2Prob(qual)/3)/used_ref_freq.freq[2]);
        case 'T': return log2((Phred2Prob(qual)/3)/used_ref_freq.freq[3]);
        default : return 0;
    }
    return 0;
}

/*
 * @abstract: calculate the pssm score.
 */ 
float pssm(const std::string & seq, const std::vector<uint8_t> & qual,
           unsigned char *alignment, int alignmentLength,
           const AGCTFreq & used_ref_freq)
{
    int seq_index=0;
    float ret_score = 0;
    bool gap_is_open=false;
    int  gap_len = 0;
    for(int i = 0 ; i<alignmentLength; i++ )
    {
        if( alignment[i] == EDLIB_EDOP_MATCH )
        {
            if(gap_is_open)
            {
                gap_is_open = false;
                ret_score -= gap_open + gap_extend*(gap_len-1);
                gap_len=0;
            }
            ret_score += base_prob_match(seq.at(seq_index),qual.at(seq_index),used_ref_freq);
            seq_index++;
        }
        else if( alignment[i] == EDLIB_EDOP_MISMATCH)
        {
            if(gap_is_open)
            {
                gap_is_open = false;
                ret_score -= gap_open + gap_extend*(gap_len-1);
                gap_len=0;
            }
            ret_score += base_prob_mismatch(seq.at(seq_index),qual.at(seq_index),used_ref_freq);
            seq_index++;
        }
        else
        {
            gap_is_open = true;
            gap_len++;
            if( alignment[i] != EDLIB_EDOP_DELETE ) seq_index++;
        }
    }
    if(gap_is_open)
    {
        gap_is_open = false;
        ret_score -= gap_open + gap_extend*(gap_len-1);
        gap_len=0;
    }
    
    return ret_score;
}

struct realign_result
{
    hts_pos_t start;
    hts_pos_t end;
    float score;
    int edit_distance;
};

/* ************************************************************************
 *
 * @section:  realign module
 *
 *************************************************************************/
//long long work = 0;
/*
 * @abstract: realign softclip reads and add info to eccs
 */
void realign_scs(std::string chrname ,ReAlignItem * item)
{
    //work ++;
    //std::cerr<<work<<std::endl;
    for( const auto & realign: item->realign_interval)
    { // foreach realign interval
        std::string forward_ref(caches.get_genome_slice(chrname,realign.next_reference_start, realign.next_reference_end));
        if ( forward_ref.size() < 1 ) continue ;
        //std::string backward_ref(seq_complementary_reverse(forward_ref));
        AGCTFreq forward_freq;// backward_freq;
        forward_freq.init(forward_ref);
        //backward_freq.init(backward_ref);
        for(const auto & s : item->scs)
        { // foreach softclip read
            if( non_colinearity(s.sc_info.headclip, s.pos_begin ,  realign.next_reference_start , realign.next_reference_end ) )
            {
                 int sc_len = s.sc_info.seqs.size();
                 int edits_allowed = adaptative_myers_k(s.sc_info.seqs.size(),edit_distance_frac);
                 EdlibAlignConfig edlib_align_config = {
                     -1,/*edits_allowed,*/ /*assumption make here*/
                     EdlibAlignMode::EDLIB_MODE_HW,
                     EdlibAlignTask::EDLIB_TASK_PATH,
                     NULL,
                     0};
                 
                 std::string used_ref;
                 AGCTFreq used_ref_freq;
                 //if( s.is_reverse )
                 //{
                 //    used_ref = backward_ref ;
                 //    used_ref_freq = backward_freq ;
                 //}
                 //else
                 //{
                     used_ref = forward_ref ;
                     used_ref_freq = forward_freq ;
                 //}
                 int hit = 0 ;
                 float min_score = int(s.sc_info.seqs.size());
                 std::vector<realign_result> r_results;
                 while( hit < nhit && min_score >= -10)
                 { // search all posible alignment
                     //std::cerr<<"call edlib 0 "<<sc_len<<std::endl;
                     EdlibAlignResult result = edlibAlign(s.sc_info.seqs.c_str(), /*query*/
                                                          sc_len,                 /*query length*/
                                                          used_ref.c_str(),       /*target*/
                                                          int(used_ref.size()),   /*target length*/
                                                          edlib_align_config);
                     //std::cerr<<"call edlib 1 "<<used_ref.size()<<std::endl;
                     if( result.status == EDLIB_STATUS_ERROR 
                       ||result.editDistance == -1 
                       ||result.numLocations < 1
                       ||result.endLocations == NULL
                       ||result.startLocations == NULL
                       ||result.alignment == NULL)
                     {
                         
                         edlibFreeAlignResult(result);
                         break; // ternimate the loop if no valid alignment
                     }
                     if( hit == 0 && result.editDistance > edits_allowed )
                     {
                         edlibFreeAlignResult(result);
                         break; // ternimate the loop if no valid alignment
                     }
                     for( int i = 0 ; i< result.numLocations ; i++ )
                     { 
                         hit++;
                         // mask aligned region
                         for( int j = result.startLocations[i]; j< result.endLocations[i]; j++ )
                             used_ref[j]='X';
                         float curr_score = pssm(s.sc_info.seqs,
                                                 s.sc_info.qual,
                                                 result.alignment,        
                                                 result.alignmentLength,  
                                                 used_ref_freq);
                         if( curr_score < min_score) min_score = curr_score;
                         //curr_score = pow(2,curr_score);
                         r_results.push_back({ result.startLocations[i] , result.endLocations[i], curr_score,result.editDistance});
                     }
                     edlibFreeAlignResult(result);
                 }
                 if( r_results.empty() || r_results[0].score <= 0 ||r_results[0].edit_distance > edits_allowed)
                     continue;
                 float prob = 1;
                 //debug
                 //if( s.qname == "V350014168L4C004R0601366836" )
                 //{
                 //    std::cerr<<prob<<std::endl;
                 //}
                 if( r_results.size() >1 )
                 {
                     double all = 0.0f;
                     for(const auto & i_result : r_results) all+=pow(2,i_result.score);
                     // we only need the best hit
                     prob = pow(2,r_results[0].score) / all;
                 }
                 //debug
                 if( prob >= prob_cutoff )
                 {
                     hts_pos_t soft_clip_start = realign.next_reference_start + r_results[0].start;
                     hts_pos_t soft_clip_end = realign.next_reference_start + r_results[0].end;
                     float score = sc_len * prob;
                     // prepare ecc info
                     eccTemp temp;
                     temp.pos_begin = s.pos_begin;
                     temp.pos_end = s.pos_end;
                     temp.leftmost_sa = soft_clip_start;
                     temp.score = round2(score);
                     temp.qname = s.qname;
                     if( s.pos_begin < soft_clip_start )
                     {
                         temp.cut_begin = s.pos_begin;
                         temp.cut_end = soft_clip_end+1;
                         item->hook_ecc(temp,true);
                     }
                     else if (s.pos_begin > soft_clip_start )
                     {
                         temp.cut_begin = soft_clip_start;
                         temp.cut_end = s.pos_end;
                         item->hook_ecc(temp,true);
                     }
                 }
            }
            else continue;
        }
    }
    return ;
}

/*
 * @abstract: realign and gather results
 */
void get_eccs_for_one_candidate_interval(std::string chrname ,ReAlignItem * item)
{
    if(! item->scs.empty() ) 
        realign_scs(chrname,item);
    if(! item->eccs.empty() )
    {
        item->assign_discordant();
    }
    else
        item->discordant_only();
}

void get_eccs()
{
    BGIQD::MultiThread::MultiThread mt;
    for(auto & itr : peaks.chromesomes)
    {
        std::string chrname = itr.first;
        region_of_chromesome & peak_array = itr.second;
        for( size_t i = 0 ; i<peak_array.size(); i++ )
        {
            ReAlignItem & item = peak_array.at(i);
            if( item.has_eccs() ) 
            {
                ReAlignItem * the_item = &item;
                mt.AddJob([chrname,the_item](){ get_eccs_for_one_candidate_interval(chrname,the_item);});
            }
        }
    }
    mt.Start(threadnum);
    //mt.Start(1);
    mt.WaitingStop();
}

void print_results(std::string output)
{
    std::ofstream ofs(output);
    //long long id = 0;
    ofs<<"chromesome\tstart\tend\tdiscordant\tsplit\tscore\n";
    for(const auto & itr : peaks.chromesomes)
    {
        const std::string & chrname = itr.first;
        const region_of_chromesome & peak_array = itr.second;
        for( size_t i = 0 ; i<peak_array.size(); i++ )
        {
            const ReAlignItem & item = peak_array.at(i);
            if( ! item.final_eccs.empty() )
            {
                //id ++;
                for(const auto & fecc : item.final_eccs)
                {
                    ofs<<chrname<<'\t'<<fecc.cut_begin<<'\t'<<fecc.cut_end<<'\t'
                             <<fecc.discordant_read_count<<'\t'<<fecc.s_qnames.size()<<'\t'
                             <<fecc.score<<'\n';
                    //if( fecc.cut_begin == 6512791&& fecc.cut_end == 6513281 ) for(const auto & x:fecc.s_qnames) std::cerr<<x<<std::endl;
                }
            }
            else if( ! item.only_discordant.empty() )
            {
                //id ++;
                for(const auto & od : item.only_discordant)
                {    ofs<<chrname<<'\t'<<od.cut_begin<<'\t'<<od.cut_end<<'\t'
                             <<od.qnames.size()<<'\t'<<0<<'\t'
                             <<0<<'\n';
                    //@debug
                    //if( od.cut_begin == 14831396 ) for(const auto & x:od.qnames) std::cerr<<x<<std::endl;
                }
            }
        }
    }
}

/* ************************************************************************
 *
 * @section: main entry.
 *
 *************************************************************************/
// prepare global datas on booting
void OnBoot()
{
    InitQual2Prob();
    hts_set_log_level(htsLogLevel::HTS_LOG_INFO);
}

void usage()
{
    std::cerr<<"\
Usage: realign_cm <-c candidate.bam>\n\
                  <-b peak.bed>\n\
                  <-g genome.fa>\n\
                  <-o output.txt>\n\
                  <-m mIS,mean value of insert size>\n\
                  <-d sIS,SD Value of insert size>\n\
                  [-S std_factor,extern realign interval by mIS+sIS*std_factor.(default 4)]\n\
                  [-t thread number(default 8)]\n\
                  [-q minimum mapping quality(default 20)]\n\
                  [-p minimum interval probability(default 0.01)]\n\
                  [-e edit distance fraction(default 0.05)]\n\
                  [-l minimum softclip length(default 8)]\n\
                  [-n nhit, maximum alignment number(default 10)]\n\
                  [-G penity for gap open(default 5)]\n\
                  [-E penity for gap extern(default 1)]\n\
                  [-P alignment probability(default 0.99)]\n";
}

int main(int argc, char ** argv)
{
    OnBoot();
    std::string bam;
    std::string peak;
    std::string genome;
    std::string output;

    int option;
    while((option = getopt(argc, argv, "c:b:g:o:m:d:S:t:q:p:e:l:n:G:E:P:h")) != -1)
    {
        switch(option)
        {
            case 'c': bam = std::string(optarg);break;
            case 'b': peak = std::string(optarg);break;
            case 'g': genome = std::string(optarg);break;
            case 'o': output = std::string(optarg);break;
            case 'm': insert_mean = atof(optarg);break;
            case 'd': insert_sd = atof(optarg);break;
            case 'S': std_fac = atoi(optarg);break;
            case 't': threadnum = atoi(optarg);break;
            case 'q': min_mapq = atoi(optarg);break;
            case 'p': interval_prob_cutoff= atof(optarg);break;
            case 'e': edit_distance_frac= atof(optarg);break;
            case 'l': min_sc_len = atoi(optarg);break;
            case 'n': nhit = atoi(optarg);break;
            case 'G': gap_open = atoi(optarg);break;
            case 'E': gap_extend = atoi(optarg);break;
            case 'P': prob_cutoff = atof(optarg);break;
            case 'h': usage();return 0;
            default:  usage();return 1;
        }
    }
    if( bam.empty() || peak.empty() || genome.empty() || output.empty() )
    {
        usage();
        return 1;
    }
    //std::cerr<<threadnum<<std::endl;
    peaks.load_peaks(peak);
    caches.load_genome(genome.c_str());
    caches.load_mates_sa_sr(bam.c_str());
    merge_candidate_mate();
    hts_log_info("merging data...");
    peaks.get_realign_intervals();
    //peaks.debug_print_realign_intervals();
    merge_candidate_sc();
    merge_candidate_split_reads();
    hts_log_info("merging data end...");
    hts_log_info("start realign softclip read, this may take a long time...");
    get_eccs();
    hts_log_info("realign softclip end...");
    print_results(output);
    hts_log_info("all done...");
    return 0;
}
