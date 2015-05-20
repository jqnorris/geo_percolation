# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:58:31 2012

@author: jqnorris
"""

import os
import multiprocessing
import subprocess
import numpy as np

run_size = 10000
number_of_runs = 13
offset = 0

arguments = range(offset, number_of_runs + offset)

for argument in arguments:
    here = os.getcwd()
    path = 'Runs/{0:d}'.format(argument)
    if not os.path.exists(path):
        os.makedirs(path)
    os.chdir(path)
    if not os.path.islink('geo_percolation'):
        subprocess.call('ln -s ../../geo_percolation', shell=True)
    os.chdir(here)

def run(argument):
    print argument
    here = os.getcwd()    
    path = 'Runs/{0:d}'.format(argument)
    os.chdir(path)
    output = 0
    output = subprocess.check_call('./geo_percolation {0}'.format(run_size),
            stdout=open(os.devnull, 'w'), shell=True)
    os.chdir(here)
    return output

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes = 1)
    results = []

    r = pool.map_async(run, arguments, callback=results.append)

    r.wait()  # Wait on the results
    
    to_file = open('job_info.txt', 'w')
    to_file.write('Job Information\n')
    to_file.write('Run Size = {0}\n'.format(run_size))
    to_file.write('Number of Runs = {0}\n'.format(number_of_runs))
    to_file.write('Offset = {0}\n'.format(offset))
    to_file.close()