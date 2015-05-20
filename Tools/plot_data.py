# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import my
import os

all_files = os.listdir('Aggregate')
data_files = [f for f in all_files if '.txt' in f]

for f_str in data_files:
    f = open('Aggregate/' + f_str, 'r')
    N = int(f.readline())
    header = f.readline()
    x = np.zeros(N, dtype=float)
    y = np.zeros(N, dtype=float)

    for i, line in enumerate(f):
        x[i], y[i] = [float(value) for value in line.split('\t')]
    
    fig, ax = my.plot_with_title()
    ax.set_title('Distribution')
    ax.plot(x, y, color=my.blue)
    my.plt.savefig('Aggregate/' + f_str.replace('.txt', '') + '.pdf', transparent=True)