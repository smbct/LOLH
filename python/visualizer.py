#!/usr/bin/python

import numpy as np
import math

import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

import histogram

# Definition of functions to produce model visual
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
def plot_rule_scores(ax, scores, selected_indexes = None, instance_info = None):

    col = ['royalblue' for _ in range(len(scores))]

    if not selected_indexes is None:
        for ind in selected_indexes:
            col[ind] = 'red'

    X = [score[0] for score in scores]
    Y = [score[1] for score in scores]

    ax.scatter(X,Y, marker='x', c=col)

    if(len(scores) > 0):
        ax.plot([scores[0][0], scores[-1][0]], [scores[0][1], scores[-1][1]], linestyle='--', linewidth=1, c='forestgreen')

    ax.set_xlabel('Positive score')
    ax.set_ylabel('Negative score')

    if instance_info is None:
        ax.set_title('Overall rule associated pareto front')
    else:
        ax.set_title('Overall rule associated pareto front for ' + instance_info[0] + '=' + str(instance_info[1]))

    return





#-------------------------------------------------------------------------------
def plot_histograms(ax, histo, normalize = True, instance_info = None):

    # normalize = True: the height of a class is relative to the total number of samples for this type of hsitogram
    # normalize = False: the height of a class is relative to the maximum number of samples in each class

    if histo.type != histogram.Histogram_type.NEGATIVE and histo.instance.n_positives() == 0:
        raise Exception('Instance has no positive samples')

    if histo.type != histogram.Histogram_type.POSITIVE and histo.instance.n_negatives() == 0:
        raise Exception('Instance has no negative samples')


    if not normalize and histo.type != histogram.Histogram_type.POSITIVE_NEGATIVE:
        normalize = True

    lim_x = 4

    error_max = len(histo.body)

    n_positives = None
    n_negatives = None

    if histo.type != histogram.Histogram_type.NEGATIVE:
        n_positives = histo.instance.n_positives()
    if histo.type != histogram.Histogram_type.POSITIVE:
        n_negatives = histo.instance.n_negatives()

    lim = ax.get_xlim()
    limy = ax.get_ylim()

    x_range = lim[1]-lim[0]
    y_range = limy[1]-limy[0]

    dx = x_range/(error_max+1)

    x_positions = [dx/2 + elt*dx for elt in range(error_max+1)]

    # plot the histogram
    for error in range(error_max+1):

        pos_height = None
        neg_height = None

        if histo.type != histogram.Histogram_type.NEGATIVE:
            pos_height = len(histo.positive_histogram[error])
            if normalize:
                pos_height/=n_positives
            else:
                pos_height/=max(n_positives, n_negatives)

        if histo.type != histogram.Histogram_type.POSITIVE:
            neg_height = len(histo.negative_histogram[error])
            if normalize:
                neg_height/=n_negatives
            else:
                neg_height/=max(n_positives, n_negatives)

        left = x_positions[error]-dx/2
        width = dx
        bottom = 0

        if histo.type != histogram.Histogram_type.NEGATIVE:
            # positive green error
            p = plt.Rectangle((left, bottom), width, pos_height, facecolor='green', alpha=0.5)
            ax.add_patch(p)

        if histo.type != histogram.Histogram_type.POSITIVE:
            # negative red error
            p = plt.Rectangle((left, bottom), width, neg_height, facecolor='red', alpha=0.5)
            ax.add_patch(p)

        if histo.type != histogram.Histogram_type.NEGATIVE:
            # positive contour
            height = pos_height*y_range
            p = plt.Rectangle((left, bottom), width, height, fill=False, color='black')
            ax.add_patch(p)

        if histo.type != histogram.Histogram_type.POSITIVE:
            # negative contour
            height = neg_height*y_range
            p = plt.Rectangle((left, bottom), width, height, fill=False, color='black')
            ax.add_patch(p)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([elt for elt in range(error_max+1)])

    ax.set_xlabel('Rule matching error on samples')
    # ax.set_xlabel('Erreur de couverture')
    ax.set_ylabel('Number of samples')
    # ax.set_ylabel('Proportion d\'états')

    # pos_patch = Patch(facecolor=colorConverter.to_rgba('green', alpha=0.5), label=r'$\mathcal{S}^+$', edgecolor='k', linewidth=0.5)
    pos_patch = Patch(facecolor=colorConverter.to_rgba('green', alpha=0.5), label='positive examples', edgecolor='k', linewidth=0.5)
    neg_patch = Patch(facecolor=colorConverter.to_rgba('red', alpha=0.5), label='negative examples', edgecolor='k', linewidth=0.5)
    ax.legend(handles=[pos_patch, neg_patch])

    ax.set_ylim((0,0.6))

    # if instance_info is None:
    #     ax.set_title('Matching score over all samples')
    # else:
    #     ax.set_title('Matching score over all samples for ' + instance_info[0] + '=' + str(instance_info[1]))

    return



#-------------------------------------------------------------------------------
def plot_global_histogram(ax, histo, instance_info = None):

    # normalize = True: the height of a class is relative to the total number of samples for this type of hsitogram
    # normalize = False: the height of a class is relative to the maximum number of samples in each class

    if histo.type != histogram.Histogram_type.GLOBAL:
        raise Exception('Can\'t display the histogram which is not global')

    lim_x = 4

    error_max = len(histo.body)

    n_samples = histo.instance.n_positives()+histo.instance.n_negatives()

    lim = ax.get_xlim()
    limy = ax.get_ylim()

    x_range = lim[1]-lim[0]
    y_range = limy[1]-limy[0]

    dx = x_range/(error_max+1)

    x_positions = [dx/2 + elt*dx for elt in range(error_max+1)]

    # plot the histogram
    for error in range(error_max+1):


        height = len(histo.positive_histogram[error])

        height /= (n_samples)

        left = x_positions[error]-dx/2
        width = dx
        bottom = 0

        # error
        p = plt.Rectangle((left, bottom), width, height, facecolor='royalblue', alpha=0.5)
        ax.add_patch(p)


        # contour
        height = height*y_range
        p = plt.Rectangle((left, bottom), width, height, fill=False, color='black')
        ax.add_patch(p)


    ax.set_xticks(x_positions)
    ax.set_xticklabels([elt for elt in range(error_max+1)])

    ax.set_xlabel('Rule matching error on samples')
    ax.set_ylabel('Number of samples')

    # if instance_info is None:
    #     ax.set_title('Matching score over all samples')
    # else:
    #     ax.set_title('Matching score over all samples for ' + instance_info[0] + '=' + str(instance_info[1]))

    return




#-------------------------------------------------------------------------------
def plot_violins(ax, df, instance, body, weights = None):

    # x positions of th violins for each gene
    violin_positions = [ind*8.5 for ind in range(len(body)) ]

    plot_index = 0
    for atom_ind in body:

        gene = instance.get_atom(atom_ind)[0]

        plot_x_pos = violin_positions[plot_index]

        positive_values = sorted(df.loc[instance._pos_samples][gene].values)
        negative_values = sorted(df.loc[instance._neg_samples][gene].values)

        ax.scatter(plot_x_pos-1.5-0.2+0.4*np.random.rand(len(positive_values)), positive_values, 1.5, marker='.', color='black', alpha=0.5)
        ax.scatter(plot_x_pos+1.5-0.2+0.4*np.random.rand(len(negative_values)), negative_values, 1.5, marker='.', color='black', alpha=0.5)

        parts = ax.violinplot([positive_values, negative_values], [plot_x_pos-1.5, plot_x_pos+1.5], points=100, widths=2., showmeans=False, showmedians=False, showextrema=False)
        ind = 0
        for pc in parts['bodies']:
            pc.set_alpha(0.5)
            pc.set_edgecolor('black')
            if ind == 0:
                pc.set_facecolor('green')
                ind += 1
            else:
                pc.set_facecolor('red')

        # parts['cbars'].set_color('black')

        plot_index += 1

    # Axes ticks and labels
    ax.set_xticks([])
    # ax.set_xticklabels([inst.get_atom(ind)[0] for ind in body])

    # Title of the plot
    ax.set_title('Solution gene value distribution across positive vs negative cells')
    ax.set_ylabel('Gene values')

    # Lengends with the colors from cell types
    legend_elements = [  Patch(facecolor=colorConverter.to_rgba('green', alpha=0.5), edgecolor='black', label='positive cells'), Patch(facecolor=colorConverter.to_rgba('red', alpha=0.5), edgecolor='black', label= 'negative cells') ]
    ax.legend(handles=legend_elements, loc='upper right')

    # Create a sub plot to show discrete normalised gene scores
    divider = make_axes_locatable(ax)
    ax_scores = divider.append_axes("bottom", 0.5, pad=0.1, sharex=ax)
    ax_scores.set_xticks(violin_positions)
    ax_scores.set_xticklabels([instance.get_atom(ind)[0] + '_' + str(instance.get_atom(ind)[1]) for ind in body])
    ax_scores.set_ylabel('Discrete scores')

    # print the gene scores

    if weights == None:
        min_score = (np.min([instance.atom_score[ind][0] for ind in range(instance.n_atoms())]), np.min([instance.atom_score[ind][1] for ind in range(instance.n_atoms())]))
        max_score = (np.max([instance.atom_score[ind][0] for ind in range(instance.n_atoms())]), np.max([instance.atom_score[ind][1] for ind in range(instance.n_atoms())]))
    else:
        gene_scores = [instance.atom_score[elt][0]*weights[0] + instance.atom_score[elt][1]*weights[1] for elt in body]
        # print('Gene scores in the solution: ', gene_scores)
        gene_scores_complete = [instance.atom_score[ind][0]*weights[0] + instance.atom_score[ind][1]*weights[1] for ind in range(instance.n_atoms())]
        min_score = np.min(gene_scores_complete)
        max_score = np.max(gene_scores_complete)

    # ax_scores.set_ylim(min_score, max_score)
    # ax_scores.ticklabel_format(axis='y', style='plain')

    for index in range(len(body)):
        atom_score = instance.atom_score[body[index]]
        pos = violin_positions[index]

        if weights == None:
            # p = plt.Rectangle((pos-1, 0.), 1, (atom_score[0]-min_score[0])/(max_score[0]-min_score[0]), facecolor=colorConverter.to_rgba('green', alpha=0.5), edgecolor='black')
            # ax_scores.add_patch(p)
            # p = plt.Rectangle((pos, 0.), 1, (atom_score[1]-min_score[1])/(max_score[1]-min_score[1]), facecolor=colorConverter.to_rgba('red', alpha=0.5), edgecolor='black')
            # ax_scores.add_patch(p)
            p = plt.Rectangle((pos-1, 0.), 1, atom_score[0]/instance.n_positives(), facecolor=colorConverter.to_rgba('green', alpha=0.5), edgecolor='black')
            ax_scores.add_patch(p)
            p = plt.Rectangle((pos, 0.), 1, atom_score[1]/instance.n_negatives(), facecolor=colorConverter.to_rgba('red', alpha=0.5), edgecolor='black')
            ax_scores.add_patch(p)
        else:
            weighted_sum = atom_score[0]*weights[0] + atom_score[1]*weights[1]
            p = plt.Rectangle((pos-1, 0.), 2, (weighted_sum-min_score)/(max_score-min_score), facecolor='royalblue', alpha=1.)
            ax_scores.add_patch(p)

    return



#-------------------------------------------------------------------------------
def create_dot_plot(ax, df, genes, cells):

    return
