#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

typedef struct
{
    unsigned short chr;
    unsigned long pos;
    signed short charge;
    float quality;
} snp;

void write_haplotype (snp* current_snp, string sample, string output_dir)
{
    stringstream output_path;
    output_path << sample << ".txt";
    fstream input_file(output_path.str(), ios::app); //make this binary ?
    if(!input_file)
    {
        cout<< "File not found, creating a new one." << endl;
    }
    input_file << current_snp->charge << endl;
    input_file.close();
}

int main () {


    return 0;
}