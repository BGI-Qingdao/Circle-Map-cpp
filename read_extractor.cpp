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
#include "htslib/thread_pool.h"
#include "htslib/hts_log.h"

/* ************************************************************************
 *
 * @section: utils functions.
 *
 *************************************************************************/
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

bool is_hardclip_read(const bam1_t *b) 
{
     uint32_t* cigar_data = bam_get_cigar(b);
     for( unsigned int i = 0 ; i< (b)->core.n_cigar; i++)
         if( bam_cigar_op(cigar_data[i]) == BAM_CHARD_CLIP )
             return true;
     return false;
}

bool bam_has_MQ(const bam1_t *b)
{
    uint8_t * data = NULL;
    data = bam_aux_get(b,"MQ");
    if ( data == NULL )
        return false;
    return true;
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
 * @abstract: ternimate the program if bam is not sorted by qname 
 */
void guarentee_bam_format(sam_hdr_t * h)
{
    kstring_t tmp = { 0, 0, NULL };
    if( sam_hdr_find_line_id(h,"HD","SO","queryname",&tmp) < 0 )
    {
        hts_log_error("input bam is not sort by qname . exit... ");
        exit(1);
    }
    free(tmp.s);
}

/* ************************************************************************
 *
 * @section: global variables.
 *
 *************************************************************************/

int mapq = 10;
int threadnum = 8;

bool NoDiscordant = false;
bool NoSoftcliped = false;
bool NoHardcliped = false;

/* ************************************************************************
 *
 * @section: bam buffer.
 *
 *************************************************************************/

struct BamBuffer
{
    samFile* the_outbam;
    bam_hdr_t * the_header;
    std::vector<bam1_t*> buffer;

    void Init(samFile* outbam, bam_hdr_t * header)
    {
        the_header = header;
        the_outbam = outbam;
    }
    void write(const bam1_t *read)
    {
        bam1_t * new_record = bam_dup1(read);
        buffer.push_back(new_record);
        //if( buffer.size() == 1024 *1024 )
        //{
        //    clean_buffer();
        //}
    }
    void clean_buffer()
    {
        for(const auto & read : buffer)
           sam_write1(the_outbam,the_header,read);
        for(bam1_t * read : buffer)
           bam_destroy1(read);
        buffer.clear();

    }
} the_buffer;


/* ************************************************************************
 *
 * @section: core functions.
 *
 *************************************************************************/

void write_cliped_read( bam1_t * read, bam1_t * mate, bool own_mapq = false)
{
    if (read->core.qual < mapq )
        return ;
    if ( ! bam_has_MQ(read) )
    {
        if( own_mapq )
            bam_aux_update_int(read,"MQ",read->core.qual);
        else
            bam_aux_update_int(read,"MQ",mate->core.qual);
    }
    if( is_softclip_read(read) && NoSoftcliped == false )
        the_buffer.write(read);
    else if ( is_hardclip_read(read) && NoHardcliped == false )
        the_buffer.write(read);
}

void LoadAndWrite(const std::string & qbam, const std::string & outbam)
{
    // try to open bam file -------------------------------------
    samFile *in = 0;
    if ((in = sam_open(qbam.c_str(), "r")) == 0)
    {
         hts_log_error("Failed to open qsorted bam : \"%s\" for reading", qbam.c_str());
         exit(1);
    }
    bam_hdr_t * header = NULL;
    header = sam_hdr_read(in);
    if (! header )
    {
        hts_log_error("Failed to load header from : \"%s\"", qbam.c_str());
        exit(1);
    }
    guarentee_bam_format(header);
    // open output bam ------------------------------------------
    samFile *out = 0;
    if ((out = sam_open(outbam.c_str(), "wb")) == 0)
    {
         hts_log_error("Failed to open outbam bam : \"%s\" for reading", outbam.c_str());
         exit(1);
    }
    // write HD:SO="unsorted" -----------------------------------
    bam_hdr_t * outheader = NULL;
    outheader = sam_hdr_dup(header);
    sam_hdr_update_hd(outheader,"SO","unsorted");
    int ret = sam_hdr_write(out,outheader);
    the_buffer.Init(out,outheader);
    // prepare thread pool --------------------------------------
    htsThreadPool p = {NULL,0};
    htsThreadPool po = {NULL,0};
    if( threadnum > 1 )
    {
        if (!(p.pool = hts_tpool_init(threadnum ) ) )
        {
            hts_log_error("Create tpool failed");
            exit(1);
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        if (!(po.pool = hts_tpool_init(threadnum ) ) )
        {
            hts_log_error("Create tpool failed");
            exit(1);
        }
        hts_set_opt(out,  HTS_OPT_THREAD_POOL, &po);
    }
    // prepare temp variables ----------------------------------
    bam1_t *bam_cache1 = bam_init1();
    bam1_t *bam_cache2 = bam_init1();
    bam1_t *curr_read1 = NULL;
    bam1_t *curr_read2 = NULL;
    bam1_t *curr_cache = bam_cache1;
    std::string r1_name;
    std::string r2_name;
    // loading bam start ---------------------------------------
    while ((sam_read1(in, header, curr_cache)) >= 0) 
    { // read one alignment from `in'
        // cache read1 and switch to next cache
        if( bam_is_r1(curr_cache) )
        {
            r1_name = std::string(bam_get_qname(curr_cache));
            curr_read1 = curr_cache ;
            if( curr_cache == bam_cache1 ) 
                curr_cache = bam_cache2;
            else
                curr_cache = bam_cache1;
            continue;
        }
        else
        {
            r2_name=std::string(bam_get_qname(curr_cache));
            if( ! bam_is_r2(curr_cache) or r1_name != r2_name)
                continue;
            curr_read2 = curr_cache;
            if( bam_is_unmap(curr_read1) or bam_is_unmap(curr_read2) )
            {// one unmap
                if( bam_is_unmap(curr_read1) == false)
                    write_cliped_read(curr_read1,curr_read2,true);
                if( bam_is_unmap(curr_read2) == false)
                    write_cliped_read(curr_read2,curr_read1,true);
            }
            else
            {
                if( bam_is_rev(curr_read2) and (!bam_is_rev(curr_read1)) )
                {
                    if(curr_read1->core.tid == curr_read2->core.tid )
                    {
                        if( curr_read2->core.pos < curr_read1->core.pos )
                        {
                            if( curr_read1->core.qual >= mapq and curr_read2->core.qual >= mapq )
                            {
                                if(NoDiscordant) continue;
                                bam_aux_update_int(curr_read1,"MQ",curr_read2->core.qual);
                                bam_aux_update_int(curr_read2,"MQ",curr_read1->core.qual);
                                the_buffer.write(curr_read1);
                                the_buffer.write(curr_read2);
                            }
                            else
                            {
                                write_cliped_read(curr_read1,curr_read2,false);
                                write_cliped_read(curr_read2,curr_read1,false);
                            }
                        }
                        else
                        { // NOT discordant
                            write_cliped_read(curr_read1,curr_read2,false);
                            write_cliped_read(curr_read2,curr_read1,false);
                        }
                    }
                    else
                    {// NOT the same chromsome
                        write_cliped_read(curr_read1,curr_read2,false);
                        write_cliped_read(curr_read2,curr_read1,false);
                    }
                }
                else
                { // NOT R2F1
                    write_cliped_read(curr_read1,curr_read2,false);
                    write_cliped_read(curr_read2,curr_read1,false);
                }
            }
        }
    }
    the_buffer.clean_buffer();
    // release resources
    bam_destroy1(bam_cache1);
    bam_destroy1(bam_cache2);
    sam_hdr_destroy(header);
    sam_hdr_destroy(outheader);
    sam_close(in);
    sam_close(out);
    if( threadnum > 1 )
    {
        hts_tpool_destroy(p.pool);
        hts_tpool_destroy(po.pool);
    }
    return ;
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
Usage: read_extractor <-i qsort.bam>\n\
                      <-o output.bam>\n\
                      [-t thread number.(Default 8)]\n\
                      [-q bwa-mem mapping quality cutoff. (Default 10)]\n\
                      [-D Turn off discordant (R2F1 oriented) read extraction. (Default not set)]\n\
                      [-S Turn off soft-clipped read extraction. (Default not set)]\n\
                      [-H Turn off hard-clipped read extraction. (Default not set)]\n";
}

int main(int argc, char ** argv)
{
    OnBoot();
    std::string bam;
    std::string output;
    int option;
    while((option = getopt(argc, argv, "i:o:t:q:hDSH")) != -1)
    {
        switch(option)
        {
            case 'h': usage(); return 0;
            case 'i': bam = std::string(optarg);break;
            case 'o': output = std::string(optarg);break;
            case 't': threadnum = atoi(optarg);break;
            case 'q': mapq = atoi(optarg);break;
            case 'D': NoDiscordant = true; break;
            case 'S': NoSoftcliped = true; break;
            case 'H': NoHardcliped = true; break;
            default : usage(); return 1;
        }
    }
    if( bam == "" or output == "" or threadnum<1 )
    {
        usage(); return 1;
    }
    LoadAndWrite(bam,output);
    return 0;
}
