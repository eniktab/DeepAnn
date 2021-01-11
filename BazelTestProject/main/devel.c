//
// Created by niktabel on 4/01/21.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define ADENINE_COEFCHR  3
#define GUANINE_COEFCHR  5

#define THYMINE_COEFCHR  7
#define CYTOSINE_COEFCHR 11

#define SNV 13

typedef struct
{
    signed short charge;
    unsigned long pos;
    float quality;
} snp;


snp* creat_haplotype(int buffer_size)
{
    snp *empty_haplotype;
    empty_haplotype = (snp *) malloc(sizeof (snp) * buffer_size);
    return empty_haplotype;
}

snp* fill_haplotype(snp* empty_haplotype, int i,
                    unsigned long c_pos, signed short c_charge, float c_quality)
{
    empty_haplotype[i].charge = c_charge;
    empty_haplotype[i].pos = c_pos;
    empty_haplotype[i].quality = c_quality;

}

signed short nuc_to_charge(char* iupac)
{
    signed short charge = 0;
    // for copy number variants
    if(iupac[0] == '<')
    {
        //<cnv[some number]
        short  cn_c;
        for (int cnt = 0; cnt<strlen(iupac); cnt++)
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
    else
    {
        for (int i = 0; i < strlen(iupac); i++)
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


signed short estimate_charge(char* reference_c, char* alternative_c)
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

int main()
{
    char* first = "A";
    char* second = "<CN0>,<CN2>,<CN3>,<CN4>";
    signed short  a =  estimate_charge(first, second);
    printf("%d", a);
    chain* A = creat_haplotype(10) ;
    //fill_haplotype(A);
    //printf("%ld \n", A->pos);
    //printf("%ld \n", A[0].pos);
    //printf("%ld \n", sizeof(ADENINE_COEFCHR));

    return 0;
};



