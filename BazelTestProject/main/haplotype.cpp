//
// Created by niktabel on 19/01/21.
//

#include "haplotype.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include <htslib/synced_bcf_reader.h>
using namespace std;

#define ADENINE_COEFCHR  3
#define GUANINE_COEFCHR  5

#define THYMINE_COEFCHR  7
#define CYTOSINE_COEFCHR 11

#define SNV 13

typedef struct
{
    unsigned short chr;
    unsigned long pos;
    signed short charge;
} snp;


signed short nuc_to_charge(const char* iupac)
{
    signed short charge = 0;

    if(iupac[0] == '<') // dealing with  copy number variants
    {
        //<cnv[some number]
        short  cn_c;
        for (long unsigned int cnt = 0; cnt<strlen(iupac); cnt++)
        {
            //printf("%c :", iupac[cnt] );
            if (cnt % 3 ==0 && iupac[cnt +1] == '>')
            {
                int cn = atoi(&iupac[cnt]);
                cn_c = cn * SNV;
                charge += cn_c;
                //printf(">> %i \n", cn);
            }
            else
            {
                int cn = atoi(&iupac[cnt]);
                cn_c = cn * SNV *10;
                charge += cn_c;
            }
        }
        return charge;
    }
    else // only snps
    {
        for (long unsigned int i = 0; i < strlen(iupac); i++)
        {
            char upper_i = toupper(iupac[i]);
            switch (upper_i)
            {
                case 'A':
                    charge += ADENINE_COEFCHR;
                    break;
                case 'G':
                    charge += GUANINE_COEFCHR;
                    break;
                case 'T':
                    charge += THYMINE_COEFCHR;
                    break;
                case 'C':
                    charge += CYTOSINE_COEFCHR;
                    break;
                default:
                    charge += SNV;
            }
        }
    }

    return charge;

}


signed short estimate_charge(const char* reference_c, const char* alternative_c)
{
    signed short charge = 0;
    if ( reference_c == alternative_c)
    {
        return charge;
    }
    signed short ref_charge = nuc_to_charge(reference_c);
    signed short alt_charge = -1 * nuc_to_charge(alternative_c);

    charge = ref_charge + alt_charge;

    return charge;

}


void write_haplotype (string output_dir, string sample, string chr, unsigned long  pos, signed short charge)
{
    stringstream output_path;
    output_path << output_dir << sample << ".txt";
    fstream input_file(output_path.str(), ios::app); //make this binary ?
    if(!input_file)
    {
        cout<< "File not found, creating a new one." << endl;
    }
    input_file << chr << ":" << pos << ":" << charge << endl;
    input_file.close();
}

void read_vcf(const char *fname)
{
    //open vcf file
    htsFile *fp  = hts_open(fname,"rb");
    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init(); //The bcf1_t structure corresponds to one VCF/BCF line
    //save for each vcf record
    while ( bcf_read(fp, hdr, rec)>=0 )
    {


        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_FMT);

        // sample names
        //TODO: is this always in the same order as the sample columns?
        //get genotype for all samples
        bcf_fmt_t *fmt = bcf_get_fmt(hdr, rec, "GT");
        if ( !fmt ) return;
        //sample genotypes
        int t =(fmt->p_len)/2;
        for (int i = 0; i<t; ++i)
        {

            //sample name
            printf(">>%s:", hdr->samples[i]);
            //read and annotate chromosome name from vcf header
            printf("%s:", bcf_hdr_id2name(hdr, rec->rid));
            //POS
            printf("%lu:", (unsigned long)rec->pos);
            //genotype
            int first_allele_p = i*2;
            int second_allele_p = first_allele_p + 1;
            int first_allele_g = bcf_gt_allele(fmt->p[first_allele_p]);
            int second_allele_g = bcf_gt_allele(fmt->p[second_allele_p]);
            signed short  charge = 0;


            //print genotype
            if (first_allele_g>0 || second_allele_g>0)
            {
                int number_alleles = (sizeof(rec->d.allele)/2);
                for (int pos = 1; pos<=number_alleles; pos++)
                {
                    charge += estimate_charge(rec->d.allele[0], rec->d.allele[pos]);
                    charge = charge/2;

                }
            }
            printf("%d \n", charge);
            write_haplotype ( "/home/niktabel/master_projects/DeepAnn/TestData/",
                              hdr->samples[i], bcf_hdr_id2name(hdr, rec->rid), (unsigned long)rec->pos, charge);
            //write_haplotype (snp* current_snp, hdr->samples[i], "/home/niktabel/master_projects/DeepAnn/TestData/";

        }

    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);


}
