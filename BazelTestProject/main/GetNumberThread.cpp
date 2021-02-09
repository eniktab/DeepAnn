//
// Created by niktabel on 10/02/21.
//

#include <omp.h>
#include <thread>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <pthread.h>
using namespace std;

int main(int argc, char** argv)
{
    stringstream output_path;
    stringstream output_path_2;

    output_path << argv[1];
    fstream out_file(output_path.str(), ios::app);

    const auto processor_count = std::thread::hardware_concurrency();
    out_file << "Number of processors: " << processor_count << endl ;
    out_file << "Number of processors available to oMP: " << omp_get_max_threads() << endl;

    return 0;
};