import matplotlib
import subprocess
import os, sys
from RunTests import run_tests

# run tests
for i in range(2, 21, 3):
    run_tests(['-o', 'O2', '-c', 'icc', '-n', str(i), '-p', '12', '-g', '32', '-m', '1', '-r', '2', '-t', '30'])
    
