# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:58:31 2012

@author: jqnorris
"""

import os
import multiprocessing
import subprocess
import numpy as np
import shutil

N = 10000000

offset = 2000
samples = 3000

arguments = range(offset, offset + samples)

for argument in arguments:
    here = os.getcwd()
    path = '{0:d}'.format(argument)
    if not os.path.exists(path):
        os.makedirs(path)
    os.chdir(path)
    if not os.path.exists('geo_percolation'):
        subprocess.call('ln -s ../geo_percolation', shell=True)
    os.chdir(here)

def run(argument):
    print argument
    here = os.getcwd()    
    path = '{0:d}'.format(argument)
    os.chdir(path)
    output = 0
    output = subprocess.check_call('./geo_percolation {0}'.format(N),
            stdout=open(os.devnull, 'w'), shell=True)
    os.chdir(here)
    return output


if __name__ == '__main__':
    pool = multiprocessing.Pool(processes = 2)
    results = []

    r = pool.map_async(run, arguments, callback=results.append)

    r.wait()  # Wait on the results
    print results

print "Done"



