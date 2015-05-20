# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 13:14:01 2013

@author: jqnorris
"""

import numpy as np
import my
import matplotlib.pyplot as plt

sim_f = open('burst_distribution.txt')
sim_N = int(sim_f.readline())
sim_header = sim_f.readline()
p_c = float(sim_header.split('=')[-1])
sim_dist = np.zeros((sim_N, 2))
for i, line in enumerate(sim_f):
    sim_dist[i] = np.fromstring(line, sep='\t')
sim_f.close()

f = open('fractures.txt', 'r')
N = int(f.readline())
header = f.readline()
data = np.zeros((N, 5))
for i, line in enumerate(f):
    data[i] = np.fromstring(line, sep="\t")
f.close()
strengths = data[:,0]

bursting = False
calc_dist = {}
for s in strengths:
    if (s < p_c):
        if bursting:
            burst_size +=1
        else:
            bursting = True
            burst_size = 1
    else:
        if burst_size in calc_dist:
            calc_dist[burst_size] += 1
        else:
            calc_dist[burst_size] = 1
        bursting = False
        burst_size = 0
        
x = sorted(calc_dist.keys())
y = [calc_dist[key] for key in x]
calc_dist = np.transpose([x,y])

fig, ax = my.plot_with_title()
ax.plot(calc_dist[:,0], calc_dist[:,1])
ax.plot(sim_dist[:,0], sim_dist[:,1])
ax.set_xscale('log')
ax.set_yscale('log')
plt.savefig('burst_check.pdf')

if(sim_dist == calc_dist).all():
    print 'Test Passed.'
else:
    print 'Test FAILED!'