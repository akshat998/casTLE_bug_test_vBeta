###############################################################################
# David Morgens
# 04/06/2016
###############################################################################
# Imports neccessary modules

from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import scipy.stats as st
import scipy.stats.mstats as ms
import csv
import sys
from collections import defaultdict
import argparse
from screenFun import *
import warnings
import math
import os


###############################################################################    
# Version number

current_version = '1.0'


###############################################################################
# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Visualizes count distribution')

# Non-optional arguments: The files containing counts, as well as an output
parser.add_argument('name', help='Name for output file.', type=str)

parser.add_argument('count_files', help='Count files.', type=str, nargs='+')

# Optional arguments:
parser.add_argument('-of', '--override_file',
                        help='Override automatic targeting to Results folder',
                        action='store_true')

parser.add_argument('-l', '--legend', type=str, nargs='+',
                        help='Name for corresponding count file')

parser.add_argument('-x', '--exclude', type=str,
                        help='Only include elements containing substrings.', nargs='+')

parser.add_argument('-s', '--search', type=str,
                        help='Use count files in indicated folder.')

parser.add_argument('-f', '--file_type', default='png',
                        help='File ending/type. Default is "png"')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input

if args.override_file:
    file_out = args.name
else:
    file_out = os.path.join('Results', args.name)

try:
    with open(file_out + '_dist.'+ args.file_type, 'w') as out_open:
        pass

    os.remove(file_out + '_dist.' + args.file_type)

except IOError:
    sys.exit('Cannot write to output file:\n' + file_out + '\n'
                + 'Use -of or --override_file to change')


###############################################################################
# Searches for other count files

if args.search:

    count_file = args.count_files[0]

    for root, dirs, files in os.walk(args.search):
        for fil in files:
            if '_counts.csv' in fil:
                args.count_files.append(os.path.join(root, fil))
                if args.legend:
                    args.legend.append(fil[: -11])


###############################################################################
# Parses count files

print('Parsing counts')

norms = []
entropies = []

# For each provided count file
for count_file in args.count_files:

    all_counts = []

    # Save the counts from each element
    with open(count_file, 'r') as count_open:

        dialect = csv.Sniffer().sniff(count_open.readline(), delimiters='\t ,')
        count_open.seek(0)
        count_csv = csv.reader(count_open, dialect)
        next(count_csv)
        for line in count_csv:

            count = False

            # Skips blank lines
            if not line or not line[0]:
                continue

            # If no exclusion characters, save line
            if not args.exclude:
                count = int(float(line[1]))

            # If exclusion character if it does not contain substring
            else:
                for ex in args.exclude:
                    if ex in line[0]:
                        count = int(line[1])
                        break

            if count:
                all_counts.append(count)


    # Normalize counts and calculates entropy
    all_counts_sorted = sorted(all_counts, reverse=True)
    tot_counts = float(sum(all_counts))

    norm_counts = [x / tot_counts for x in all_counts_sorted]
    norms.append(norm_counts)

    entropy = 0

    for count in all_counts:
        entropy += -1 * (count / tot_counts) * math.log(count / tot_counts, 2)

    norm_entropy = entropy / math.log(len(all_counts),2)
    entropies.append(entropy)

# Finds the total number of elements and normalizes entropies to this value
tot_elements = max(map(len, norms))
norm_entropies = [ent / math.log(tot_elements, 2) for ent in entropies]

print('Total elements: ' + str(tot_elements))


###############################################################################
# Plots the figure

print('Plotting figure')

fig = plt.figure(dpi=400)
ax = plt.subplot(111)

# Plots all distributions on single plot
for norm in norms:
    ax.plot(range(len(norm)), norm)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

# Adds legend
if args.legend:

    fontP = FontProperties()
    fontP.set_size('small')
    new_legend = [name + '; diversity = ' + str(sigDig(ent)) \
                        for name, ent in zip(args.legend, norm_entropies)]
    plt.legend(new_legend, prop=fontP, loc='best', bbox_to_anchor=(1, 1))

# Plots in log scale
plt.yscale('log')
plt.xlabel('Elements')
plt.ylabel('Normalized Counts')
plt.title('Distribution of counts')

plt.savefig(file_out + '_dist.' + args.file_type)
plt.close()


###############################################################################
