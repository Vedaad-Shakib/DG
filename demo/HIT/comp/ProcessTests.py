import sys, getopt # to parse the arguments
import subprocess # to run bash commands
import os
import re
import time
import csv
import filecmp
import datetime

TEST_DIRECTORY = "/home/vshakib/LYDG/demo/HIT/comp"
BIN_DIRECTORY = "/home/vshakib/LYDG/build/bin"

def main(argv):
    opts, args = getopt.getopt(argv, "j:o:c:n:p:g:m:t:")

    if len(opts) < 8:
        print "all 8 parameters are required"
        print "usage: python ProcessTests.py -o O2 -c gcc -n 5 -p 12 -g 8 -m 1 -r 1 -t 10"
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
        else:
            print "%s flag not recognized" % (arg)

    print "job %s has finished" % (job_number)

    elapsed, cpu = extract_time(job_number)
    print "elapsed time: %s" % (elapsed)
    print "cpu time: %s" % (cpu)

    today = datetime.datetime.today().strftime("%d/%m")
    write_csv(job_number, "HIT", today, str(n_nodes), str(n_processors), str(n_timesteps), str(grid_size),
              str(compiler), str(optimization), str(mv2_enable_affinity), elapsed, cpu)
    write_profile(job_number)
    

def write_profile(job_number):
    subprocess.call("cd %s; gprof %s/plydg gmon.out > profiles/profile.%s" % (TEST_DIRECTORY, BIN_DIRECTORY, job_number)

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

if __name__ == "__main__":
    main(sys.argv[1:])
