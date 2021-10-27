#!/usr/bin/python

import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.ticker as ticker
from matplotlib import cm

import sys
sys.path.append('../../python')

from instance import Instance

import visualizer
import histogram

# 0: gene name
# 1: gene discrete value
# 2: quality
# 3 : number of positive examples
# 4 : number of negative examples

#-------------------------------------------------------------------------------
def read_quality_file(filename):

    file = open(filename, 'r')
    content = file.read().splitlines()
    file.close()

    data = []

    for line in content:
        tokens = line.split(' ')
        data.append((tokens[0], int(tokens[1]), float(tokens[2]), int(tokens[3]), int(tokens[4])))

    return data

#-------------------------------------------------------------------------------
def plot_quality(data, ax):

    ax.set_xlabel('Quality')
    ax.set_ylabel(r'Proportion between $\mathcal{S}^+$ and $\mathcal{S}^-$')

    indexes = [ind for ind in range(len(data)) if data[ind][2] > 0]

    qualities = [data[ind][2] for ind in indexes]
    proportions = [data[ind][3]/(data[ind][3]+data[ind][4]) for ind in indexes]

    colors = ['red' if data[ind][1] == 0 else 'blue' for ind in indexes]

    ax.scatter(qualities, proportions, s=2, marker='o', c=colors)

    # ind2 = 0
    # for ind in indexes:
    #     if qualities[ind2]  >= 0.8 and math.fabs(proportions[ind2]-0.5) <= 0.4:
    #         ax.text(qualities[ind2], proportions[ind2], data[ind][0] + '_' + str(data[ind][1])).set_clip_on(True)
    #     ind2 += 1

    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0]) / (ax.get_ylim()[1]-ax.get_ylim()[0]))

    return

print('hello dynamics')

# filename = 'coexpression_quality.txt'
# # filename = 'transitions_quality.txt'
# data = read_quality_file(filename)
# fig, ax = plt.subplots()
# plot_quality(data, ax)
# plt.show()

filename = '../../dataset/Imagine/transitions.csv'
df = pd.read_csv(filename, index_col=0)
print(df.head())

successors = {}

for index in df.index:
    left = df['T-1'][index]
    right = df['T'][index]

    if left in successors:
        successors[left].append(right)
    else:
        successors[left] = [right]

print(successors)

nsuc = [len(successors[elt]) for elt in successors]
nsuc.sort()
print(nsuc)
