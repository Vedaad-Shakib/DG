import sys, getopt # to parse the arguments
import subprocess # to run bash commands
import os
import re
import time
import csv
import filecmp
import datetime

BASE_DIRECTORY = "/home/vshakib/LYDG"
TEST_DIRECTORY = BASE_DIRECTORY+"/demo/HIT/comp"
BIN_DIRECTORY = BASE_DIRECTORY+"/build/bin"
BUMP_DIRECTORY = BASE_DIRECTORY+"/demo/HIT/comp/bump"
SCRIPT_DIRECTORY = BASE_DIRECTORY+"/demo/HIT/comp/lydg"
PROFILE_DIRECTORY = BASE_DIRECTORY+"/demo/HIT/comp/profile"
STDERR_DIRECTORY = BASE_DIRECTORY+"/demo/HIT/comp/stderr"
STDOUT_DIRECTORY = BASE_DIRECTORY+"/demo/HIT/comp/stdout"

def main(argv):
    opts, args = getopt.getopt(argv, "j:o:c:n:p:g:m:t:i:")

    if len(opts) < 9:
        print "all 9 parameters are required"
        print "usage: python ProcessTests.py -j job_id -o O2 -c gcc -n 5 -p 12 -g 8 -m 1 -r 1 -t 10 -i <tmp_id>"
        sys.exit(2)

    for opt, arg in opts:
        if arg == "":
            print "%s requires an option to be set" % (opt)
            sys.exit(2)
        elif opt=="-j":
            job_number = re.findall("\d+", arg)[0]
        elif opt=="-o":
            optimization = arg
        elif opt=="-c":
            compiler = arg
        elif opt=="-n":
            n_nodes = arg
        elif opt=="-p":
            n_processors = arg
        elif opt=="-g":
            grid_size = arg
        elif opt=="-m":
            mv2_enable_affinity = arg
        elif opt=="-t":
            n_timesteps = arg
        elif opt=="-i":
            tmp_id = arg
        else:
            print "%s flag not recognized" % (arg)

    print "job %s has finished" % (job_number)

    elapsed, cpu = extract_time(job_number)
    print "elapsed time: %s" % (elapsed)
    print "cpu time: %s" % (cpu)

    today = datetime.datetime.today().strftime("%m/%d")
    write_csv(job_number, "HIT", today, str(n_nodes), str(n_processors), str(n_timesteps), str(grid_size),
              str(compiler), str(optimization), str(mv2_enable_affinity), elapsed, cpu, str(n_timesteps))
    #write_profile(job_number)

    rename_tmp_files(job_number, tmp_id)

# renames tmp job and script files with the job id    
def rename_tmp_files(job_number, tmp_id):
    subprocess.call("mv %s/bump.%s.job %s/bump.%s.job" % (BUMP_DIRECTORY, tmp_id,
                                                          BUMP_DIRECTORY, job_number), shell=True)
    subprocess.call("mv %s/lydg.%s.script %s/lydg.%s.script" % (SCRIPT_DIRECTORY, tmp_id,
                                                                SCRIPT_DIRECTORY, job_number), shell=True)
    subprocess.call("mv %s/%s.certainty-fe.stanford.edu %s/stderr.%s" % (STDERR_DIRECTORY, job_number,
                                                                         STDERR_DIRECTORY, job_number), shell=True)
    subprocess.call("mv %s/%s.certainty-fe.stanford.edu %s/stdout.%s" % (STDOUT_DIRECTORY, job_number,
                                                                         STDOUT_DIRECTORY, job_number), shell=True)

# generates a profile using gprof and writes it to a directory
# not used because gprof is kind of useless    
def write_profile(job_number):
    subprocess.call("cd %s; gprof %s/plydg gmon.out > %s/profile.%s" % (TEST_DIRECTORY, BIN_DIRECTORY, PROFILE_DIRECTORY, job_number))

def write_csv(*arg):
    fout = open("baseline.csv", "ab")
    fout.write(",".join(arg))
    fout.write("\n")
    fout.close()

# extracts cpu_time and elapsed_time from log file    
def extract_time(job_number):
    current_dir = os.path.dirname(os.path.abspath(__file__)) 
    log_name = current_dir+"/logs/%s.certainty-fe.stanford.edu" % (job_number)
    log = open(log_name, "r").read()

    log = log.split("\n")

    elapsed_time = find_start(log, "Elapsed time")
    elapsed_time = re.findall("\d+.\d+", elapsed_time)[0]
    
    cpu_time = find_start(log, "CPU time")
    cpu_time = re.findall("\d+.\d+", cpu_time)[0]

    return elapsed_time, cpu_time

# returns the first element in array "arr" that starts with "start"
def find_start(arr, start):
    for i in arr:
        if i.startswith(start):
            return i
    return None

if __name__ == "__main__":
    main(sys.argv[1:])
