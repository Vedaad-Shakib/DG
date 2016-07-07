### usage: python runTests.py -o O2 -c gcc -n 5 -p 12 -g 8 -m 1 -r 3 -t 10
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
BUMP_DIRECTORY = "/home/vshakib/LYDG/demo/HIT/comp/bump"
SCRIPT_DIRECTORY = "/home/vshakib/LYDG/demo/HIT/comp/lydg"

def main(argv):
    script_template=open("lydg_template.script", "r").read()
    bump_template=open("bump_template.job", "r").read()
    cmake_template=open("%s/CMakeLists_template.txt" % (CMAKE_DIRECTORY), "r").read()
    prev_settings=open("%s/PrevSettings.txt" % (CMAKE_DIRECTORY), "r").read().split("\n")

    opts, args = getopt.getopt(argv, "o:c:n:p:g:m:r:t:")
    
    job_numbers = []
    job_status = []

    n_jobs_finished = 0

    devnull = open(os.devnull, "w")

    rebuild = False # whether to rebuild the code or not


    if len(opts) < 8:
        print "all 8 parameters are required"
        print "usage: python RunTests.py -o O2 -c icc -n 10 -p 12 -g 32 -m 1 -r 2 -t 10"
        sys.exit(2)

    print "initialized script"
        
    for opt, arg in opts:
        if arg == "":
            print opt+" requires an option to be set"
            sys.exit(2)
        elif opt=="-o":
            print "optimization: %s" % (arg)

            opt_dict = {"O2": "RelWithDebInfo", "O3": "Release"}
            if prev_settings[0] != arg:
                rebuild = True
                cmake_template = cmake_template.replace("OPTIMIZATION_LEVEL", opt_dict[arg])
                
            script_template = script_template.replace("OPTIMIZATION", arg)
            optimization = arg
        elif opt=="-c":
            print "compiler: %s" % (arg)
            
            if prev_settings[1] != arg:
                rebuild = True
                subprocess.call("export CC=$(which %s)" % (arg), shell=True)
                subprocess.call("Export CXX=$(which %s)" % (arg), shell=True)
            
            script_template = script_template.replace("COMPILER", arg)
            compiler = arg
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

            script_template = script_template.replace("GRID_SIZE", arg)
            bump_template = bump_template.replace("GRID_SIZE", arg)
            grid_size = arg
        elif opt=="-m":
            print "mv2_enable_affinity: %s" % (arg)
            
            script_template = script_template.replace("MV2_SETTING", arg)
            mv2_enable_affinity = arg
        elif opt=="-r":
            print "# runs: %s" % (arg)

            script_template = script_template.replace("N_RUNS", arg)
            n_runs = int(arg)
        elif opt=="-t":
            print "$ timesteps: %s" % (arg)

            script_template = script_template.replace("N_TIMESTEPS", arg)
            bump_template = bump_template.replace("N_TIMESTEPS", arg)
            n_timesteps = arg
        else:
            print "%s flag not recognized" % (arg)

    script_template = script_template.replace("TOT_PROCESSORS", str(n_processors*n_nodes))

    tmp_id = str(hash(str(opts))) # generate unique identifier for each set of options; will be renamed once job id is assigned
    bump_run=open("%s/bump.%s.job" % (BUMP_DIRECTORY, tmp_id), "w")
    bump_run.write(bump_template)
    bump_run.close()

    script_template = script_template.replace("BUMP_FILE", "%s/bump.%s.job" % (BUMP_DIRECTORY, tmp_id))
    script_template = script_template.replace("TMP_ID", tmp_id)
    script_run=open("%s/lydg.%s.script" % (SCRIPT_DIRECTORY, tmp_id), "w")
    script_run.write(script_template+"\n")
    script_run.close()

    cmake_run=open("%s/CMakeLists.txt" % (CMAKE_DIRECTORY), "w")
    cmake_run.write(cmake_template)
    cmake_run.close()

    curr_settings=open("%s/PrevSettings.txt" % (CMAKE_DIRECTORY), "w")
    curr_settings.write(optimization+"\n")
    curr_settings.write(compiler)
    curr_settings.close()

    print "customized job file"

    #if rebuild:
        #make_project()

    '''
    # run jobs
    for i in range(n_runs):
        job = subprocess.check_output("qsub %s/lydg.%s.script" % (SCRIPT_DIRECTORY, tmp_id), shell=True)
        job_numbers.append(str(re.findall("^\d*", str(job))[0]))
        job_status.append(-1)
        print "running job %s" % (str(job_numbers[i]))
    '''

    print "job numbers: %s" % (", ".join(job_numbers))

    print "all jobs have been submitted "                

def make_project():    
    print "rebuilding with adjusted optimization settings"
    subprocess.call("cd %s; rm -rf *; cmake ../; make install" % (BUILD_DIRECTORY), shell=True)

# returns the first element in array "arr" that starts with "start"
def find_start(arr, start):
    for i in arr:
        if i.startswith(start):
            return i
    return None

if __name__ == "__main__":
   main(sys.argv[1:])

