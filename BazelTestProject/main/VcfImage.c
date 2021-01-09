//compile:
//$gcc -c -o test-vcf.o test-vcf.c
//$gcc -pthread   -o test-vcf test-vcf.o ../htslib/libhts.a -lz
//should change ../htslib/libhts.a to your own path to libhts.a
//https://github.com/samtools/htslib/blob/develop/htslib/vcf.h#L222

#include <stdio.h>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include <htslib/synced_bcf_reader.h>


void read_vcf(char *fname)
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


            if (first_allele_g>0)
            {
                printf("%s", rec->d.allele[1]);
            }
            else
            {
                printf("%s", rec->d.allele[0]);
            }
            if (second_allele_g>0)
            {
                printf("%s\n", rec->d.allele[1]);
            }
            else
            {
                printf("%s\n", rec->d.allele[0]);
            }



        }


/*


        //rec->d.allele[0] is REF other is ALT
        printf("%i\n", *rec->d.allele[0]);

        for (int i=1; i<rec->n_allele; ++i)
        {
            printf("%s\n", rec->d.allele[i]);
        }

 */

    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);


}

int main()
{
    read_vcf("/home/niktabel/master_projects/DeepAnn/TestData/head-trim");
    return 0;
};
