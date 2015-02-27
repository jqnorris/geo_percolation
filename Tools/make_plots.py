# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:58:31 2012

@author: jqnorris
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import plots
plt.rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)

N = 10000000

samples = 5000

arguments = range(samples)

main_path = '3D Combined'

branching_distribution = {}
burst_distribution = {}
mass_l_distribution = {}
mass_r_distribution = {}

# If run data has not be aggregated
if not os.path.exists(main_path):
    
    # For each run
    for argument in arguments:
        
        # Change directory to run folder
        here = os.getcwd()
        path = '{0:d}'.format(argument)
        os.chdir(path)
        
        # Read in branching statistics
        f = open('branch_distribution.txt', 'r')
        max_order = f.readline()
        header = f.readline()        
        for line in f:
            i, j, length, counts = [int(i) for i in line.split('\t')]
            if (i, j, length) in branching_distribution:              
                branching_distribution[(i, j, length)] += counts
            else:
                branching_distribution[(i, j, length)] = counts
        f.close()
        
        # Read in burst statistics
        f = open('burst_distribution.txt', 'r')
        num_sizes = f.readline()
        header = f.readline()
        for line in f:
            size, counts = [int(i) for i in line.split('\t')]
            if size in burst_distribution:
                burst_distribution[size] += counts
            else:
                burst_distribution[size] = counts
        f.close()
        
        # Read in chemical level statistics
        f = open('mass_l_distribution.txt', 'r')
        num_levels = f.readline()
        header = f.readline()
        for line in f:
            level, counts = [int(i) for i in line.split('\t')]
            if level in mass_l_distribution:
                mass_l_distribution[level] += counts
            else:
                mass_l_distribution[level] = counts
        f.close()
        
        # Read in mass statistics
        f = open('mass_r_distribution.txt', 'r')
        max_radius = f.readline()
        header = f.readline()
        for line in f:
            radius, counts = [int(i) for i in line.split('\t')]
            if radius in mass_r_distribution:
                mass_r_distribution[radius] += counts
            else:
                mass_r_distribution[radius] = counts
        f.close()
        
        os.chdir(here)
    
    # Make folder for agregated statistics
    here = os.getcwd()
    os.makedirs(main_path)
    os.chdir(main_path)
        
    # Write branching statistics to file
    f = open('branch_distribution.txt', 'w')
    max_order = sorted(branching_distribution)[-1][0]
    f.write('{}\n'.format(max_order))
    f.write('Aggregate Branching Statistics\n')
    for key in sorted(branching_distribution):
        f.write('{0}\t{1}\t{2}\t{3}\n'.format(key[0], key[1], key[2], branching_distribution[key]))
    f.close()
    
    # Write burst statistics to file
    f = open('burst_distribution.txt', 'w')
    f.write('{}\n'.format(len(burst_distribution)))
    f.write('Aggregate Burst Distribution\n')
    for key in sorted(burst_distribution):
        f.write('{0}\t{1}\n'.format(key, burst_distribution[key]))
    f.close()
    
    # Write chemical level statistics to file
    f = open('mass_l_distribution.txt', 'w')
    f.write('{}\n'.format(len(mass_l_distribution)))
    f.write('Aggregate Chemical Level Statistics\n')
    for key in sorted(mass_l_distribution):
        f.write('{0}\t{1}\n'.format(key, mass_l_distribution[key]))
    f.close()
    
    # Write mass statistics to file
    f = open('mass_r_distribution.txt', 'w')
    f.write('{}\n'.format(len(mass_r_distribution)))
    f.write('Aggregate Mass Statistics\n')
    for key in sorted(mass_r_distribution):
        f.write('{0}\t{1}\n'.format(key, mass_r_distribution[key]))
    f.close()
    
    os.chdir(here)
else:
    # Load data from file
    here = os.getcwd()
    os.chdir(main_path)
        
    # Read in branching statistics
    f = open('branch_distribution.txt', 'r')
    max_order = f.readline()
    header = f.readline()        
    for line in f:
        i, j, length, counts = [int(i) for i in line.split('\t')]
        if (i, j, length) in branching_distribution:              
            branching_distribution[(i, j, length)] += counts
        else:
            branching_distribution[(i, j, length)] = counts
    f.close()
        
    # Read in burst statistics
    f = open('burst_distribution.txt', 'r')
    num_sizes = f.readline()
    header = f.readline()
    for line in f:
        size, counts = [int(i) for i in line.split('\t')]
        if size in burst_distribution:
            burst_distribution[size] += counts
        else:
            burst_distribution[size] = counts
    f.close()
        
    # Read in chemical level statistics
    f = open('mass_l_distribution.txt', 'r')
    num_levels = f.readline()
    header = f.readline()
    for line in f:
        level, counts = [int(i) for i in line.split('\t')]
        if level in mass_l_distribution:
            mass_l_distribution[level] += counts
        else:
            mass_l_distribution[level] = counts
    f.close()
    
    # Read in mass statistics
    f = open('mass_r_distribution.txt', 'r')
    max_radius = f.readline()
    header = f.readline()
    for line in f:
        radius, counts = [int(i) for i in line.split('\t')]
        if radius in mass_r_distribution:
            mass_r_distribution[radius] += counts
        else:
            mass_r_distribution[radius] = counts
    f.close()
        
    os.chdir(here)

# Make plots in aggregate directory
here = os.getcwd()
os.chdir(main_path)

# Make plot of number-order and length-order branching
horton_strahler_distribution = {}

for key in branching_distribution:
    new_key = (key[0], key[2])
    if new_key in horton_strahler_distribution:
        horton_strahler_distribution[new_key] += branching_distribution[key]
    else:
        horton_strahler_distribution[new_key] = branching_distribution[key]

number_order = {}
        
for key in horton_strahler_distribution:
    if key[0] in number_order:
        number_order[key[0]] += horton_strahler_distribution[key]
    else:
        number_order[key[0]] = horton_strahler_distribution[key]

number_order = np.array(sorted(number_order.items()))

x = number_order[:,0]
y = number_order[:,1]/float(samples)
plt.figure(figsize=(5,4))
plt.plot(x, y)
plt.xlabel('Branch Order')
plt.ylabel('Number')
plt.yscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('number_order.pdf')

length_order = {}

for key in horton_strahler_distribution:
    if key[0] in length_order:
        length_order[key[0]].append((key[1],  horton_strahler_distribution[key]))
    else:
        length_order[key[0]] = [(key[1],  horton_strahler_distribution[key]),]

for key in length_order:
    data = np.array(length_order[key])
    length_order[key] = np.average(data[:,0], weights=data[:,1])
    
length_order = np.array(sorted(length_order.items()))
x = length_order[:,0]
y = length_order[:,1]
plt.figure(figsize=(5,4))
plt.plot(x, y)
plt.xlabel('Branch Order')
plt.ylabel('Average Length')
plt.yscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('length_order.pdf')

# Make plot of Tokunaga statistics
Tokunaga_distribution = {}

for key in branching_distribution:
    new_key = (key[0], key[1])
    if new_key in Tokunaga_distribution:
        Tokunaga_distribution[new_key] += branching_distribution[key]
    else:
        Tokunaga_distribution[new_key] = branching_distribution[key]

order = sorted(Tokunaga_distribution)[-1][0]
Tokunaga = np.zeros((order, order))

for key in Tokunaga_distribution:
    Tokunaga[key[0]-1][key[1]-1] = Tokunaga_distribution[key]/float(number_order[key[1]-1][1])

k_array = np.arange(1, order)
T_k_array = 0.0*k_array

for k in k_array:
    for i in range(0, order-k):
        T_k_array[k-1] += Tokunaga[i][i+k]
    T_k_array[k-1] += T_k_array[k-1]/float(order-k)

plt.figure(figsize=(5,4))
plt.plot(k_array, T_k_array)
plt.xlabel('k')
plt.ylabel('T_k')
plt.yscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('tokunaga_branching.pdf', trasparent=True)   

# Make plot of burst distribution
burst_distribution = np.array(sorted(burst_distribution.items()), dtype=np.int64)
x = burst_distribution[:,0]
y = burst_distribution[:,1]
quick_binned = y[:-1]/(1.0*(x[1:]-x[:-1]))
plt.figure(figsize=(5,4))
plt.plot(x[:-1], quick_binned)
plt.xlabel('Burst Size')
plt.ylabel('Counts')
plt.yscale('log')
plt.xscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('burst_distribution.pdf', trasparent=True)

# Make plot of chemical level masses
mass_l_distribution = np.array(sorted(mass_l_distribution.items()), dtype=np.int64)
x = mass_l_distribution[:,0]
y = mass_l_distribution[:,1]
sum_y = np.cumsum(y)/float(samples)
plt.figure(figsize=(5,4))
plt.plot(x, sum_y)
plt.xlabel('Chemical Level')
plt.ylabel('Cummulative Mass')
plt.xscale('log')
#plt.yscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('mass_l_distribution.pdf', trasparent=True)

# Make plot of mass as a function of radius
mass_r_distribution = np.array(sorted(mass_r_distribution.items()), dtype=np.int64)
x = mass_r_distribution[:,0]
y = mass_r_distribution[:,1]
plots.plot_distribution([x, y/(float(samples)*float(N))])
sum_y = np.cumsum(y)/float(samples)
plt.figure(figsize=(5,4))
plt.plot(x, sum_y)
plt.xlabel('Radius')
plt.ylabel('Cummulative Mass')
#plt.xscale('log')
#plt.yscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('mass_r_distribution.pdf', trasparent=True)

# Make plot of the local slopes of the mass as a function of radius
plt.figure(figsize=(5,4))
plt.plot(x[1:], sum_y[1:] - sum_y[:-1])
plt.xlabel('Radius')
plt.ylabel('Cummulative Mass')
#plt.xscale('log')
#plt.yscale('log')
plt.subplots_adjust(bottom=0.2, left=0.2)
plt.savefig('diff_mass_r_distribution.pdf', trasparent=True)

# Plot Branch length distributions
##plt.figure()
##
##cmap = plt.get_cmap('jet')
##
##current_order = 0
##
##for key in sorted(branching_distribution):
##    if key[0] != current_order:
##        plt.savefig('test_{0}.pdf'.format(current_order), transparent=True)
##        plt.figure()
##        plt.yscale('log')
##        current_order = key[0]
##    plt.plot([key[2]], [branching_distribution[key]], marker='.', markeredgecolor='none', markerfacecolor=cmap(key[1]/float(order)))  


os.chdir(here)

