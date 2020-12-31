#include <stdio.h>
#include <iostream> 
#include <fstream>
#include <vector>
#include <string>
//#include "read_sam.h"
#include "cram.h"


using namespace std;

int main()
 {
	samFile *fp_in = hts_open("/home/niktabel/master_projects/SampleGenomicData/ref.cram", "r");
	bam_hdr_t *cram_hdr = sam_hdr_read(fp_in);
	hts_idx_t *cram_indx = sam_index_load(fp_in, "/home/niktabel/master_projects/SampleGenomicData/ref.cram.crai");
	bam1_t *aln = bam_init1(); //initialize an alignment
	hts_itr_t *iter = sam_itr_querys(cram_indx, cram_hdr, "chm1:10-15");
	
	//if (iter == NULL) 
	//{ // region invalid or reference name not found
	//	printf(stderr, "region specifies an invalid region or unknown reference. Continue anyway.\n");
	//}
	//else
	{
		while ((sam_itr_next(fp_in, iter, aln)) >= 0) {

			int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
			char *chr = cram_hdr->target_name[aln->core.tid] ; //contig name (chromosome)
			uint32_t len = aln->core.l_qseq; //length of the read.
			
			uint8_t *q = bam_get_seq(aln); //quality string
			uint32_t q2 = aln->core.qual ; //mapping quality
			
			char *qseq = (char *)malloc(len);
	
			for(int i=0; i< len ; i++)
			{
				qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
			}
	
			printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);
			}	
	}
	
	bam_destroy1(aln);
	hts_idx_destroy(cram_indx);
	sam_close(fp_in);
    
	return 0;
 };
