// compile:  g++ ProtoExample.cpp Genome.pb.cc `pkg-config --cflags --libs protobuf` -o protoexample
#include <fstream>
#include <iostream>
#include <string>
#include "Genome.pb.h"
using namespace std;


void store_variant(Genome::Variant* variant)
{
    cout << "chr?";
    short int chr;
    cin >> chr;
    variant->set_chromosome(chr);

    cout << "pos?";
    long pos;
    cin >> pos;
    variant->set_pos(pos);

    cout << "charge?";
    short int charge;
    cin >> charge;
    variant->set_charge(charge);
}

int main () {

    Genome::Person haplotype;

    fstream input("file.txt", ios::in | ios::binary);
    if (!input) {
        cout << "File not found, creating a new one." << endl;
    } else if (!haplotype.ParseFromIstream(&input)) {
        cout << "Failed to open file." << endl;
        return -1;
    }

    store_variant(haplotype.add_haplotype());
    fstream output("file.txt", ios::out | ios::binary);
    if (!haplotype.SerializeToOstream(&output)) {
        cerr << "Failed to write the variant record" << endl;
        return -1;
    }

    //Delete all global objects allocated by libprotobuf.
    google::protobuf::ShutdownProtobufLibrary();


    return 0;
}