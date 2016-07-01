### usage: python runTests.py -o O2 -c gcc -n 5 -p 12 -g box8 -m 1 -r 3 -t 10
#
# -o optimization                 O2, O3
# -c compiler                     gcc, icc
# -n number of nodes              >= 1
# -p number of processors         1-12
# -g grid size                    8, 16, 32, 64
# -m mv2_enable_affinity          0, 1
# -r number of runs               > 1
# -t number of timesteps          > 1

import sys, getopt # to parse the arguments
import subprocess # to run bash commands
import os
import re
import time
import csv
import filecmp
import datetime

CMAKE_DIRECTORY = "/home/vshakib/LYDG"
BUILD_DIRECTORY = "/home/vshakib/LYDG/build"

def main(argv):
    script_template=open("lydg_template.script", "r").read()
    script_run=open("lydg_run.script", "w")

    bump_template=open("bump_template.job", "r").read()
    bump_run=open("bump_run.job", "w")
    
    opts, args = getopt.getopt(argv, "o:c:n:p:g:m:r:t:")

    job_numbers = []
    job_status = []

    n_jobs_finished = 0

    devnull = open(os.devnull, "w")

    rebuild = False # whether to rebuild the code or not

    print "initialized script"

    if len(opts) < 8:
        print "all 8 parameters are required"
        print "usage: python RunTests.py -o O2 -c gcc -n 5 -p 12 -g 8 -m 1 -r 1 -t 10"
        sys.exit(2)
        
    for opt, arg in opts:
        if arg == "":
            print opt+" requires an option to be set"
            sys.exit(2)
        elif opt=="-o":
            print "optimization: %s" % (arg)
            
            if arg=="O3":
                if not filecmp.cmp('%s/CMakeLists_Release.txt' % (CMAKE_DIRECTORY),
                                   '%s/CMakeLists.txt' % (CMAKE_DIRECTORY)):
                    subprocess.call("cp %s/CMakeLists_Release.txt %s/CMakeLists.txt" %
                                    (CMAKE_DIRECTORY, CMAKE_DIRECTORY), shell=True)
                    rebuild = True
            elif arg=="O2":
                if not filecmp.cmp('%s/CMakeLists_RelWithDebInfo.txt' % (CMAKE_DIRECTORY),
                                   '%s/CMakeLists.txt' % (CMAKE_DIRECTORY)):
                    subprocess.call("cp %s/CMakeLists_RelWithDebInfo.txt %s/CMakeLists.txt" %
                                    (CMAKE_DIRECTORY, CMAKE_DIRECTORY), shell=True)
                    rebuild = True
            else:
                print "invalid argument %s for option %s" % (arg, opt)
                
            if rebuild:
                print "rebuilding with adjusted optimization settings"
                subprocess.call("cd %s; rm -rf *; cmake ../; make install" % (BUILD_DIRECTORY), shell=True)
            
            optimization = arg
        elif opt=="-c":
            # not functional
            print "compiler: %s" % (arg)
            
            compiler = arg
            pass
        elif opt=="-n":
            print "# nodes: %s" % (arg)
            script_template = script_template.replace("N_NODES", arg)
            
            n_nodes = int(arg)
        elif opt=="-p":
            print "# processors per node: %s" % (arg)
            script_template = script_template.replace("N_PROCESSORS", arg)
            
            n_processors = int(arg)
        elif opt=="-g":
            print "grid size: %s" % (arg)
            bump_template = bump_template.replace("GRID_SIZE", arg)
            
            grid_size = arg
        elif opt=="-m":
            print "mv2_enable_affinity: %s" % (arg)
            script_template = script_template.replace("MV2_SETTING", arg)

            mv2_enable_affinity = arg
        elif opt=="-r":
            print "# runs: %s" % (arg)
            
            n_runs = int(arg)
        elif opt=="-t":
            print "$ timesteps: %s" % (arg)
            bump_template = bump_template.replace("N_TIMESTEPS", arg)

            n_timesteps = arg
        else:
            print "%s flag not recognized" % (arg)

    script_template = script_template.replace("TOT_PROCESSORS", str(n_processors*n_nodes))
    
    script_run.write(script_template+"\n")
    script_run.close()

    bump_run.write(bump_template)
    bump_run.close()

    print "customized job"
    
    # run jobs
    for i in range(n_runs):
        job = subprocess.check_output("qsub lydg_run.script", shell=True)
        job_numbers.append(str(re.findall("^\d*", str(job))[0]))
        job_status.append(-1)
        print "running job %s" % (str(job_numbers[i]))

    print "job numbers: %s" % (", ".join(job_numbers))

    time.sleep(5)
    
    # check to see if jobs are finished
    # job status: -1-queued, 0-running, 1-finished
    while n_jobs_finished < n_runs:
        time.sleep(5)
        for i in range(n_runs):
            if job_status[i] != 1: # if queued or running, check if done
                exit_status = subprocess.call("checkjob %s" % (job_numbers[i]),
                                              shell=True, stdout=devnull, stderr=devnull)
                if exit_status == 1:
                    print "job %s has finished" % (job_numbers[i])
                    job_status[i] = 1
                    n_jobs_finished += 1

                    elapsed, cpu = extract_time(job_numbers[i])
                    print "elapsed time: %s" % (elapsed)
                    print "cpu time: %s" % (cpu)

                    today = datetime.datetime.today().strftime("%d/%m")
                    write_csv(job_numbers[i], "HIT", today, str(n_nodes), str(n_processors), str(n_timesteps), str(grid_size),
                              str(compiler), str(optimization), str(mv2_enable_affinity), elapsed, cpu)
            if job_status[i] == -1: # if still queued, check if running
                check_job = subprocess.check_output("checkjob %s" % (job_numbers[i]), shell=True)
                if "State: Running" in check_job:
                    print "job %s has queued" % (job_numbers[i])
                    job_status[i] = 0
                
    print "\nall jobs have finished"                

def write_csv(*arg):
    fout = open("baseline.csv", "ab")
    fout.write(",".join(arg))
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

