#!/usr/bin/python

import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import os


#-------------------------------------------------------------------------------
def read_quality_file(filename):

    # 0: gene name
    # 1: gene discrete value
    # 2: quality
    # 3 : number of positive examples
    # 4 : number of negative examples

    file = open(filename, 'r')
    content = file.read().splitlines()
    file.close()

    data = []

    for line in content:
        tokens = line.split(' ')
        data.append((tokens[0], int(tokens[1]), float(tokens[2]), int(tokens[3]), int(tokens[4])))

    return data

#-------------------------------------------------------------------------------
def plot_quality(data, ax, annotations):

    indexes = [ind for ind in range(len(data)) if data[ind][2] > 0]

    if not annotations:
        qualities = [data[ind][2] for ind in indexes]
        proportions = [data[ind][3]/(data[ind][3]+data[ind][4]) for ind in indexes]
        colors = ['firebrick' if data[ind][1] == 0 else 'royalblue' for ind in indexes]
        ax.scatter(qualities, proportions, s=0.6, marker='.', c=colors)
        # ax.scatter(qualities, proportions, s=6, marker='.', c=colors)
    else:
        qualities_pos = [data[ind][2] for ind in indexes if data[ind][1] > 0]
        proportions_pos = [data[ind][3]/(data[ind][3]+data[ind][4]) for ind in indexes if data[ind][1] > 0]
        qualities_nul = [data[ind][2] for ind in indexes if data[ind][1] == 0]
        proportions_nul = [data[ind][3]/(data[ind][3]+data[ind][4]) for ind in indexes if data[ind][1] == 0]
        ax.scatter(qualities_nul, proportions_nul, s=6, marker='.', c='firebrick', label='null atoms')
        ax.scatter(qualities_pos, proportions_pos, s=6, marker='.', c='royalblue', label='positive atoms')


        ax.set_xlabel('correlations quality')


        ax.set_ylabel(r'proportion between $\mathcal{S}^+$ and $\mathcal{S}^-$')
        ax.legend(loc='lower left', markerscale=4.0)

    # ind2 = 0
    # for ind in indexes:
    #     if qualities[ind2]  >= 0.8 and math.fabs(proportions[ind2]-0.5) <= 0.4:
    #         ax.text(qualities[ind2], proportions[ind2], data[ind][0] + '_' + str(data[ind][1])).set_clip_on(True)
    #     ind2 += 1


    ax.set_xlim((0.5, 1.))
    ax.set_ylim((0., 1.))

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
def plot_tau(subfig, rho, delta, files_parameters, config = 0, labels = 0):

    path = '../../dataset/Imagine/transitions_parameters_expe'

    if config == 0:
        nrows = 3
        ncols = 2
    else:
        nrows = 5
        ncols = 1

    axs = subfig.subplots(nrows, ncols)
    axs = axs.flat

    if config == 0:
        subfig.delaxes(axs[-1])

    subfig.suptitle('$\\rho = ' + str(rho) + ', \\delta = ' + str(delta) + '$' )

    selected_param = [elt for elt in files_parameters if elt[1][1] == rho and elt[1][2] == delta]
    selected_param.sort(key=lambda elt: elt[1][0])

    index = 0
    for elt in selected_param:

        data = read_quality_file(path + '/' + elt[0])
        plot_quality(data, axs[index], False)

        if config == 0:
            if elt[1][0] < 0:
                title = '$\\tau = \epsilon$'
            else:
                title = '$\\tau = ' + str(elt[1][0]) + '$'

            axs[index].set_title(title)

        axs[index].xaxis.set_major_locator(ticker.NullLocator())
        axs[index].yaxis.set_major_locator(ticker.NullLocator())

        index += 1

    return


#-------------------------------------------------------------------------------
def plot_transitions_parameters_quality():

    # list all files in a folder

    path = '../../dataset/Imagine/transitions_parameters_expe'
    files = os.listdir(path)

    files_parameters = extract_files_parameters(files)

    fig = plt.figure(constrained_layout=True, figsize=(10, 4))
    fig.suptitle('Evaluation of the transition parameters values (part 1)')

    subfigs = fig.subfigures(1, 2, wspace=0.15)

    rho = 0
    delta = 1
    plot_tau(subfigs[0], rho, delta, files_parameters, 0, 0)

    rho = 1
    delta = 1
    plot_tau(subfigs[1], rho, delta, files_parameters, 0, 0)

    #####################################################################
    fig = plt.figure(constrained_layout=True, figsize=(10, 4))

    fig.suptitle('Evaluation of the transition parameters values (part 2)')

    # subfigs = fig.subfigures(1, 4, wspace=0.15)
    subfigs = fig.subfigures(1, 5, wspace=0.05)

    subfigs[0].suptitle(' ')
    axs = subfigs[0].subplots(5, 1)
    axs = axs.flat

    ind = 0
    for tau in ['\\epsilon', '0.2', '0.4', '0.6', '0.8']:
        axs[ind].text(0.4, 0.5, '$\\tau = ' + tau + '$', size='large', horizontalalignment='left', verticalalignment='center')
        axs[ind].axis('off')
        ind += 1

    rho = 2
    delta = 1
    plot_tau(subfigs[1], rho, delta, files_parameters, 1, 1)

    rho = 0
    delta = 2
    plot_tau(subfigs[2], rho, delta, files_parameters, 1, 0)

    rho = 1
    delta = 2
    plot_tau(subfigs[3], rho, delta, files_parameters, 1, 0)

    rho = 2
    delta = 2
    plot_tau(subfigs[4], rho, delta, files_parameters, 1, 0)

    plt.show()

    return

#-------------------------------------------------------------------------------
def selected_parameters_quality():

    filename = 'transitions_quality.txt'
    data = read_quality_file(filename)
    fig, ax = plt.subplots()

    ax.set_title(r'Corelation quality with $(\tau = 0.6, \rho = 0, \delta = 1)$')
    plot_quality(data, ax, True)


    plt.show()

    return

# plot_transitions_parameters_quality()

# filename = 'coexpression_quality.txt'
# data = read_quality_file(filename)
# fig, ax = plt.subplots()
#
# plot_quality(data, ax, True)
#
# plt.show()
