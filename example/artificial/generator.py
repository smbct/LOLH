#!/usr/bin/python

import math

import numpy as np
import pandas as pd

# random generators
import random
from numpy.random import default_rng

from scipy.stats import t # student distribution

# matplotlib
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
def print_ascii_histograms(selected_atoms, atoms, positives, negatives):

    pos_histogram = [0 for _ in range(len(selected_atoms)+1)]
    neg_histogram = [0 for _ in range(len(selected_atoms)+1)]

    for pos in positives:
        score = 0
        for ind_atom in selected_atoms:
            atom = atoms[ind_atom]
            if pos[atom[0]] != atom[1]:
                score += 1
        pos_histogram[score] += 1

    for neg in negatives:
        score = 0
        for ind_atom in selected_atoms:
            atom = atoms[ind_atom]
            if neg[atom[0]] != atom[1]:
                score += 1
        neg_histogram[score] += 1


    for score in range(len(pos_histogram)):
        pos_histogram[score] /= float(len(positives))

    for score in range(len(neg_histogram)):
        neg_histogram[score] /= float(len(negatives))

    # print('positive histogram: ')
    # print(pos_histogram)
    #
    # print('negative histogram:')
    # print(neg_histogram)

    height = 10
    output = ''

    for y in range(height, 0, -1):
        for score in pos_histogram:
            if score*height*2 >= y:
                output += '*'
            else:
                output += ' '
        output += '\n'
    for score in pos_histogram:
        output += '.'
    print(output)

    output = ''
    for y in range(height, 0, -1):
        for score in neg_histogram:
            if score*height*2 >= y:
                output += '*'
            else:
                output += ' '
        output += '\n'
    for score in neg_histogram:
        output += '.'
    print(output)

    return

#-------------------------------------------------------------------------------
def compute_atom_scores(nVar, sample_score, sample_uniform, nPos, nNeg):

    a = nNeg
    b = -nPos

    atoms = [(ind_var, 0) for ind_var in range(nVar)]
    atoms += [(ind_var, 1) for ind_var in range(nVar)]

    atom_scores = []

    for ind_atom in range(len(sample_score)):
        c = sample_score[ind_atom]*nPos*nNeg
        if c < 0:
            p1 = (-c/nPos, 0)
            p2 = (nPos, nPos+c/nNeg)
        else:
            p1 = (0, c/nPos)
            p2 = (nNeg-c/nPos, nNeg)
        v = (p2[0]-p1[0], p2[1]-p1[1])
        lateral = sample_uniform[ind_atom]
        p_star = (p1[0]+v[0]*lateral, p1[1]+v[1]*lateral)
        p_star = (math.floor(p_star[0]), math.floor(p_star[1]))
        atom_scores.append(p_star)

    return atoms, atom_scores

#-------------------------------------------------------------------------------
def generate_positives_negatives(atom_scores, nVar, nPos, nNeg):

    positives = [ [0 for _ in range(nVar)] for _ in range(nPos) ]
    negatives = [ [0 for _ in range(nVar)] for _ in range(nNeg) ]

    for ind_atom in range(nVar):
        score = atom_scores[ind_atom]

        # positive examples
        pos_indexes = [ind for ind in range(nPos)]
        np.random.shuffle(pos_indexes)
        for ind_pos in range(nPos):
            if ind_pos < score[0]:
                positives[pos_indexes[ind_pos]][ind_atom] = 1
            else:
                positives[pos_indexes[ind_pos]][ind_atom] = 0

        # negative examples
        neg_indexes = [ind for ind in range(nNeg)]
        np.random.shuffle(neg_indexes)
        for ind_neg in range(nNeg):
            if ind_neg < score[1]:
                negatives[neg_indexes[ind_neg]][ind_atom] = 1
            else:
                negatives[neg_indexes[ind_neg]][ind_atom] = 0

    return positives, negatives

#-------------------------------------------------------------------------------
def generate_dataset(nVariables, nPositives, nNegatives, filename):

    nPos = 1000
    nNeg = 1000

    nVar = 10000
    nAtoms = nVar*2

    random_seed = 46
    # random_seed = 52 # worst dataset

    random.seed(random_seed)
    np.random.seed(random_seed)

    # sampling over the normal distribution
    rng = default_rng(random_seed)

    # mu, sigma = 0, 0.1 # mean and standard deviation
    # sample_normal = rng.normal(mu, sigma, nAtoms)

    sample_uniform = rng.uniform(0, 1, nVar*2)

    # generate the atom scores between 1 and -1 (generate only the score of one of the two atoms)
    # (the other one is symmetrical)
    df = 10
    sample_score = t.rvs(df, size=nVar)

    # normal is normalized
    max_normal = np.max(sample_score)
    min_normal = np.min(sample_score)

    # print(min_normal, max_normal)

    max_abs = max(max_normal, -min_normal) + 0.2
    for ind in range(len(sample_score)):
        sample_score[ind] = sample_score[ind]/max_abs
        sample_score[ind] *= 0.80


    # plot the distributions of the generated values
    # fig, axs = plt.subplots(2)
    # axs[0].hist(sample_score, edgecolor='k', bins=100)
    # axs[1].hist(sample_uniform, edgecolor='k', bins=100)


    # generate the atoms and the positive and negative errors of the atoms based on their scores
    atoms, atom_scores = compute_atom_scores(nVar, sample_score, sample_uniform, nPos, nNeg)

    # generating positive and negative examples
    positives, negatives = generate_positives_negatives(atom_scores, nVar, nPos, nNeg)

    # create the symmetrical atoms
    atom_scores_sym = []
    for score in atom_scores:
        atom_scores_sym.append( (nPos-score[0], nNeg-score[1]) )
    atom_scores += atom_scores_sym

    # plot the atom errors
    fig,ax = plt.subplots()
    ax.scatter([elt[0] for elt in atom_scores], [elt[1] for elt in atom_scores], marker='x')
    ax.plot([0, nPos], [0, nNeg], color='k', label='score=0.0')
    ax.set_xlim([0, nPos])
    ax.set_ylim([0, nNeg])
    ax.legend(loc='lower right')
    ax.set_xlabel('atom positive error')
    ax.set_ylabel('atom negative error')
    ax.set_title('positive and negative errors for all the atoms')
    ax.set_aspect('equal')


    plt.show()

    # compute the best 20 atoms to create a rule
    ind_atoms_sorted = [ind for ind in range(len(atom_scores))]
    score = [-nNeg*error[0] +nPos*error[1] for error in atom_scores]
    ind_atoms_sorted.sort(key = lambda ind: score[ind], reverse=True)
    selected_atoms = ind_atoms_sorted[:20]

    # plot the scores of the atoms as colors
    # fig, ax = plt.subplots()
    # ax.scatter([atom_scores[ind][0] for ind in ind_atoms_sorted], [atom_scores[ind][1] for ind in ind_atoms_sorted], c=[score[ind] for ind in ind_atoms_sorted], marker='x')
    # ax.set_xlim([0, nPos])
    # ax.set_ylim([0, nNeg])
    # ax.plot([0, nPos], [0, nNeg], color='red')

    # plot histograms of the rule with the best 20 atoms
    # print_ascii_histograms(selected_atoms, atoms, positives, negatives)

    # create a csv file with the matrix (ordered positive and negative examples)
    matrix = positives + negatives
    variables = ['v_'+str(ind) for ind in range(nVar)]
    examples = ['s_' + str(ind) for ind in range(nPos+nNeg)]

    # exportation into a dataframe
    # global dataframe
    dataframe = pd.DataFrame(matrix, columns = variables, index = examples)
    dataframe.to_csv(filename)



    return

# generate_dataset(10000, 1000, 1000, '')
