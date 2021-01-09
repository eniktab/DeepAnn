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

signed short estimate_charge(char reference_c[], char alternative_c[])
{
    int charge = 0;
    for (int i = 0; i < strlen(reference_c); i++)
    {
        char upper_i = toupper(reference_c[i]);
        switch (upper_i) {
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
        }
    }
        for (int i = 0; i < strlen(alternative_c); i++)
        {
            char upper_i = toupper(alternative_c[i]);
            switch(upper_i)
            {
                case 'A':
                    charge -= ADENINE_COEFCHR;
                    break;
                case 'G':
                    charge -= GUANINE_COEFCHR;
                    break;
                case 'T':
                    charge -= THYMINE_COEFCHR;
                    break;
                case 'C':
                    charge -= CYTOSINE_COEFCHR ;
                    break;
            }
        }
    return charge;

}

int main()
{
    int X = 1;
    signed short  a =  estimate_charge("att", "atC");
    printf("%d", a);
    //snp* A = creat_haplotype(10) ;
    //fill_haplotype(A);
    //printf("%ld \n", A->pos);
    //printf("%ld \n", A[0].pos);
    //printf("%ld \n", sizeof(ADENINE_COEFCHR));

    return 0;
};



