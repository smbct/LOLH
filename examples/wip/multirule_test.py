#!/usr/bin/python

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors as mcolors

import sys
sys.path.append('../../python')

from instance import Instance
from solver import Solver
import visualizer
import histogram

#------------------------------------------------------------------------------
# main function
#------------------------------------------------------------------------------
def celltype_classification():

    print('\n')
    print('Comparison of LOLH and PRIDE algorithms')
    print('\n')

    # read the discrete matrix
    filename = '../../dataset/Imagine/discrete_matrix.csv'
    df_discrete = pd.read_csv(filename, index_col=0)
    # print(df_discrete.head())

    # read the seurat cell types
    df_celltypes = pd.read_csv('../../dataset/Imagine/cell_types.csv', index_col = 0)
    df_celltypes.rename(columns={'cellType_final': 'Label'}, inplace=True)
    # print(df_celltypes.head())

    # initialization of the classification instance: classification of the NK cells
    celltype = 'NK'
    instance = Instance.create_cluster_instance(df_discrete.copy(deep=False), df_celltypes, celltype)

    print('Classification of the NK cells')
    print('- ', instance.n_positives(), ' positive examples')
    print('- ', instance.n_negatives(), ' negative examples')
    print('\n')

    solver = Solver(instance)

    # create a rule from a threshold
    sorted_atoms, scores = solver.select_best_atoms_fast()
    threshold = 0.55
    print('threshold: ', threshold)

    # print(scores[:10])
    threshold *= instance.n_positives()*instance.n_negatives()
    lolh_body = [sorted_atoms[ind] for ind in range(len(sorted_atoms)) if scores[ind] >= threshold] # 21 atoms selected

    print('-> Output of LOLH: ', [instance.get_atom(index)[0] + '_' + str(instance.get_atom(index)[1]) for index in lolh_body])


    print('Visualisations of the LOLH rule')

    fig, ax = plt.subplots()
    ax.scatter([instance.atom_score[ind][0] for ind in range(instance.n_atoms())], [instance.atom_score[ind][1] for ind in range(instance.n_atoms())], marker='x', alpha=0.5)
    ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], c='k')

    ax.set_xlim([0, instance.n_positives()])
    ax.set_ylim([0, instance.n_negatives()])

    ax.set_xlabel('positive error')
    ax.set_ylabel('negative error')

    ax.set_title('NK cell classification')

    ax.set_aspect(instance.n_positives()/instance.n_negatives())

    threshold = 0.5

    # score = threshold
    A = (0, threshold*instance.n_negatives())
    B = (instance.n_positives()*(1.-threshold), instance.n_negatives())

    print(A, B)

    ax.plot([A[0], B[0]], [A[1], B[1]], '--', color='red')

    print('split the positive examples into two sets of ', instance.n_positives()/2, ' examples')

    # compute the new worst score possible by removing all atoms in half

    plt.show()

celltype_classification()
