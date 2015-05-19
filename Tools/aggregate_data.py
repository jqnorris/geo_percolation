import numpy as np
import os

# Helper functions
def number(some_str):
    if '.' in some_str:
        return float(some_str)
    else:
        return int(some_str)

params = {}
# Load job information
f = open('job_info.txt', 'r')
header = f.readline()
for line in f:
    key, str_value = line.split('=')
    params[key.strip()] = number(str_value)

# Compare job information to file structure
from_file_structure = sorted(map(int, os.listdir('Runs')))
if 'Number of Runs' in params:
    number_of_runs = params['Number of Runs']
    offset = 0
    if 'Offset' in params:
        offset = params['Offset']
    from_job_information = range(offset, offset + number_of_runs)
else:
    print 'Warning: Job information is incomplete.'
    print 'Using file structure.'
    run_id_list = from_file_structure
    number_of_runs = len(run_id_list)  
if from_file_structure != from_job_information:
    print 'Warning: Job information does not match the file structure.'
    if set(from_job_information).issubset(from_file_structure):
        print 'Using job information.'
        run_id_list = from_job_information
    else:
        print 'WARNING! Job information and file structure incompatible.'
        quit()
else:
    print 'File structure consistent with job information.'
    run_id_list = from_job_information
        
# Determine data files to aggregate
all_files = os.listdir('Runs/' + str(run_id_list[0]))
data_files = [f_str for f_str in all_files if 'distribution' in f_str]

agg_data = {}
for f_str in data_files:
        agg_data[f_str] = {}
data_headers = {}
# Load data from files and aggregate
for i in run_id_list:
    here = os.getcwd()
    os.chdir('Runs/' + str(i) + '/')
    for f_str in data_files:
        if f_str in os.listdir(os.getcwd()):
            f = open(f_str, 'r')
        else:
            print 'Warning: Missing data file.'
            continue
        number_of_values = f.readline()
        header = f.readline()
        if f_str in data_headers:
            if header != data_headers[f_str]:
                print 'Warning: Inconsistent data files.'
        else:
            data_headers[f_str] = header
        for line in f:
            key, value_str = line.split('\t')
            if key in agg_data[f_str]:
                agg_data[f_str][key] += number(value_str)
            else:
                agg_data[f_str][key] = number(value_str)
        f.close()
    os.chdir(here)

# Write aggregate data to file
here = os.getcwd()
if not os.path.exists('Aggregate'):
    os.makedirs('Aggregate')
os.chdir('Aggregate')
for f_str in data_files:
    f = open(f_str, 'w')
    f.write(str(len(agg_data[f_str])) + '\n')
    f.write(data_headers[f_str])
    sorted_keys = sorted([int(key) for key in agg_data[f_str].keys()])
    for key in sorted_keys:
        f.write(str(key) + '\t' + str(agg_data[f_str][str(key)]/float(number_of_runs)) + '\n')
    f.close()

    
