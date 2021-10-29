#!/usr/bin/python

import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.ticker as ticker
from matplotlib import cm

import os

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

    # ax.set_xlabel('Quality')
    # ax.set_ylabel(r'Proportion between $\mathcal{S}^+$ and $\mathcal{S}^-$')

    indexes = [ind for ind in range(len(data)) if data[ind][2] > 0]

    qualities = [data[ind][2] for ind in indexes]
    proportions = [data[ind][3]/(data[ind][3]+data[ind][4]) for ind in indexes]

    colors = ['firebrick' if data[ind][1] == 0 else 'royalblue' for ind in indexes]

    ax.scatter(qualities, proportions, s=0.6, marker='.', c=colors)

    ax.set_xlim((0.5, 1.))
    ax.set_ylim((0., 1.))

    # ind2 = 0
    # for ind in indexes:
    #     if qualities[ind2]  >= 0.8 and math.fabs(proportions[ind2]-0.5) <= 0.4:
    #         ax.text(qualities[ind2], proportions[ind2], data[ind][0] + '_' + str(data[ind][1])).set_clip_on(True)
    #     ind2 += 1

    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0]) / (ax.get_ylim()[1]-ax.get_ylim()[0]))

    return

#-------------------------------------------------------------------------------
def extract_files_parameters(files):

    # select files from the experiment (.txt)
    expe_files = [file for file in files if file[-3:] == 'txt']

    result = []

    # extract parameters
    for file in expe_files:
        tokens = file[:-4].split('_')
        parameters = (float(tokens[3]), int(tokens[5]), int(tokens[7]))
        result.append((file, parameters))

    return result

#-------------------------------------------------------------------------------
def plot_tau(subfig, rho, delta, files_parameters):

    nrows = 3
    ncols = 2

    axs = subfig.subplots(nrows, ncols)
    axs = axs.flat
    # fig.tight_layout()
    subfig.delaxes(axs[-1])

    subfig.suptitle('$\\rho = ' + str(rho) + ', \\delta = ' + str(delta) + '$' )

    selected_param = [elt for elt in files_parameters if elt[1][1] == rho and elt[1][2] == delta]
    selected_param.sort(key=lambda elt: elt[1][0])
    # print(selected_param)

    index = 0
    for elt in selected_param:

        data = read_quality_file(path + '/' + elt[0])
        plot_quality(data, axs[index])

        # title = '$\\tau = ' + str(elt[1][0]) + ', \\rho = ' + str(elt[1][1]) + ', \\delta = ' + str(elt[1][2]) + '$'
        title = '$\\tau = ' + str(elt[1][0]) + '$'

        # print(title)
        axs[index].set_title(title)

        axs[index].xaxis.set_major_locator(ticker.NullLocator())
        axs[index].yaxis.set_major_locator(ticker.NullLocator())

        index += 1

    return


# list all files in a folder

path = '../../dataset/Imagine/transitions_parameters_expe'
files = os.listdir(path)

files_parameters = extract_files_parameters(files)



fig = plt.figure(constrained_layout=True, figsize=(10, 4))
subfigs = fig.subfigures(1, 2, wspace=0.15)

rho = 0
delta = 1
plot_tau(subfigs[0], rho, delta, files_parameters)

rho = 1
delta = 1
plot_tau(subfigs[1], rho, delta, files_parameters)

#################################################################

fig = plt.figure(constrained_layout=True, figsize=(10, 4))
subfigs = fig.subfigures(1, 2, wspace=0.15)

rho = 2
delta = 1
plot_tau(subfigs[0], rho, delta, files_parameters)

rho = 0
delta = 2
plot_tau(subfigs[1], rho, delta, files_parameters)

#####################################################################

fig = plt.figure(constrained_layout=True, figsize=(10, 4))
subfigs = fig.subfigures(1, 2, wspace=0.15)

rho = 1
delta = 2
plot_tau(subfigs[0], rho, delta, files_parameters)

rho = 2
delta = 2
plot_tau(subfigs[1], rho, delta, files_parameters)


plt.show()
