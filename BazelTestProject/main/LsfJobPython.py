import subprocess
import sys
import os
import time
import cProfile


def prepare_io(list_of_files, exe_file, input_path, output_path, job_number):
    # read file names
    with open(list_of_files, "r") as files_to_read:
        list_files = files_to_read.read().split("\n")
    job_number = int(job_number) - 1
    input_file = list_files[job_number]
    output_dir = os.path.join(output_path, input_file).replace(".vcf.gz", "/")
    zip_output_path = os.path.join(output_path, input_file).replace(".vcf.gz", ".tar.xz")
    to_read = os.path.join(input_path, input_file)

    if not os.path.isdir(output_dir):
        subprocess.run("mkdir {}".format(output_dir), shell=True, stdout=subprocess.PIPE)
    with open(os.path.join(output_path, input_file).replace(".vcf.gz", "_logs.txt"), "a") as log_file:
        log_file.write("{} \n".format(input_file))
        # run vcf to tensor -- c++ code
        exe = "{} {} {}".format(exe_file, to_read, output_dir)
        start = time.time()
        subprocess.run(exe, shell=True, stdout=subprocess.PIPE)
        end = time.time()
        elapsed = (end - start) / 360
        log_file.write("{} was done in {} hours \n".format(exe, elapsed))

        # zip output files
        exe_2 = "tar -cjf {} {}".format(zip_output_path, output_dir)
        start = time.time()
        subprocess.run(exe_2, shell=True, stdout=subprocess.PIPE)
        end = time.time()
        elapsed = (end - start) / 360
        log_file.write("{} was done in {} hours \n".format(exe_2, elapsed))

        # remove residual files

        exe_3 = "rm -r {}".format(output_dir)
        log_file.write("{} started \n".format(exe_3))
        subprocess.run(exe_3, shell=True, stdout=subprocess.PIPE)
        log_file.write("{} was done \n".format(exe_3))


def main(argv):
    if int(argv[4]) == 1:
        prof = cProfile.Profile()
        prof.enable()
        prepare_io(list_of_files=argv[0], exe_file=argv[1], input_path=argv[2], output_path=argv[3], job_number=argv[4])
        prof.disable()
        prof.dump_stats(os.path.join(argv[2], "sample.prof"))
    else:
        prepare_io(list_of_files=argv[0], exe_file=argv[1], input_path=argv[2], output_path=argv[3], job_number=argv[4])


if __name__ == "__main__":
    main(sys.argv[1:])
