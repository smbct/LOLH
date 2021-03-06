#!/usr/bin/python3

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from matplotlib import cm
from matplotlib import colors as mcolors

# LOLH objects

import sys
sys.path.append('../../python')

from instance import Instance
from solver import Solver

import visualizer
import histogram

#-------------------------------------------------------------------------------
def NK_classification():

    # import the cell types
    file_name = '../../dataset/Imagine/cell_types.csv'
    df_cell_types = pd.read_csv(file_name, index_col=0)
    df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)

    # import the UMAP 2d representation
    file_name = '../../dataset/Imagine/umap_coordinates.csv'
    umap_coordinates = pd.read_csv(file_name, index_col = 0)

    # import the single cell discretized matrix
    file_name = '../../dataset/Imagine/discrete_matrix.csv'
    df = pd.read_csv(file_name, index_col=0)

    # creation of the classification instance
    cell_type = 'NK'
    instance = Instance.create_cluster_instance(df.copy(deep=False), df_cell_types, cell_type)

    positive_cells = instance._pos_samples
    negative_cells = instance._neg_samples

    score_thresold = 0.55

    # compute the body from the score threshold
    solver = Solver(instance)
    LOLH_rule, _ = solver.select_best_atoms_threshold(score_thresold)


    ################################################################################
    # display the atom scores only
    # fig, ax = plt.subplots()
    # ax.scatter([instance.atom_score[atom_index][0] for atom_index in range(instance.n_atoms())], [instance.atom_score[atom_index][1] for atom_index in range(instance.n_atoms())], marker='x', alpha=0.7)
    # ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='black', label='score=0.0', zorder=3)
    # ax.set_xlabel('erreur positive')
    # ax.set_ylabel('erreur négative')
    # ax.legend(loc='lower right')
    # ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))



    ################################################################################
    # display the atom scores
    fig, ax = plt.subplots()

    other_atoms_indexes = [atom_index for atom_index in range(instance.n_atoms()) if not atom_index in LOLH_rule]

    # ax.scatter([instance.atom_score[atom_index][0] for atom_index in LOLH_rule], [instance.atom_score[atom_index][1] for atom_index in LOLH_rule], marker='x', color='red', label='atomes sélectionnés', zorder=2)

    ax.scatter([instance.atom_score[atom_index][0] for atom_index in LOLH_rule], [instance.atom_score[atom_index][1] for atom_index in LOLH_rule], marker='x', color='red', label='selected atoms', zorder=2)


    ax.scatter([instance.atom_score[atom_index][0] for atom_index in other_atoms_indexes], [instance.atom_score[atom_index][1] for atom_index in other_atoms_indexes], marker='x', zorder=1, label='other atoms', alpha=0.7)

    # ax.scatter([instance.atom_score[atom_index][0] for atom_index in other_atoms_indexes], [instance.atom_score[atom_index][1] for atom_index in other_atoms_indexes], marker='x', zorder=1, label='autres atomes', alpha=0.7)

    # equation score = 0
    ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='black', label='score=0.0', zorder=3)

    # equation score = score_threshold (0.55)
    A = (0,score_thresold*instance.n_negatives())
    B = ((1.-score_thresold)*instance.n_positives(),instance.n_negatives())
    ax.plot([A[0], B[0]], [A[1], B[1]], '--', color='black', label='score='+str(score_thresold), zorder=3)

    ax.set_aspect(instance.n_positives()/instance.n_negatives())

    ax.set_xlim((0, instance.n_positives()))
    ax.set_ylim((0, instance.n_negatives()))
    # ax.set_xlabel('erreur positive')
    ax.set_xlabel('positive error')
    # ax.set_ylabel('erreur négative')
    ax.set_ylabel('negative error')

    ax.set_title('Atom errors for the classification of the NK cells by LOLH')

    ax.legend(loc='lower right')





    ################################################################################
    # display the histograms

    fig,ax = plt.subplots()
    histograms = histogram.Histogram(instance, LOLH_rule)
    visualizer.plot_histograms(ax, histograms, True)
    ax.set_title('Histograms of the matching error for the classification of NK cells')
    ax.set_ylim((0, 0.3))




    ################################################################################
    # display the UMAP with LOLH rule score on the cells

    cell_scores = {index:instance.compute_rule_error(LOLH_rule, index) for index in umap_coordinates.index}
    sorted_cell_indexes = list(umap_coordinates.index)
    sorted_cell_indexes.sort(key=lambda index: cell_scores[index], reverse=True)
    umap_coordinates = umap_coordinates.loc[sorted_cell_indexes]

    # fig, ax = plt.subplots(dpi=400)
    fig, ax = plt.subplots()

    plasma = cm.get_cmap('plasma_r', 256)
    cnorm = mcolors.Normalize(vmin=0, vmax=len(LOLH_rule)+1)

    ax.scatter(umap_coordinates['UMAP_1'][:], umap_coordinates['UMAP_2'][:], cmap=plasma, norm=cnorm, c=[cell_scores[index] for index in umap_coordinates.index], s=2)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_aspect( (ax.get_xlim()[1]-ax.get_xlim()[0]) / (ax.get_ylim()[1]-ax.get_ylim()[0])  )

    ax.set_title('Matching error of the cells with LOLH rule (NK classification)')

    cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
    cbar.set_label('Matching error')

    umap_limits = (ax.get_xlim(), ax.get_ylim())





    ################################################################################
    # selection of the cells that are not NK, with a matching error <= 15
    # negative_cells_selected = [elt for error in range(16) for elt in histograms.negative_histogram[error]]
    # negative_cells_selected.reverse()
    #
    # other_cells = [index for index in umap_coordinates.index if not index in negative_cells_selected and df_cell_types['Label'][index] != 'NK']
    #
    # cnorm = mcolors.Normalize(vmin=0, vmax=15)
    #
    # fig,ax = plt.subplots()
    # ax.set_xlabel('UMAP 1')
    # ax.set_ylabel('UMAP 2')
    # ax.scatter(umap_coordinates.loc[other_cells]['UMAP_1'][:], umap_coordinates.loc[other_cells]['UMAP_2'][:], c='grey', s=2, alpha=0.2)
    # ax.scatter(umap_coordinates.loc[negative_cells_selected]['UMAP_1'][:], umap_coordinates.loc[negative_cells_selected]['UMAP_2'][:], c=[cell_scores[cell] for cell in negative_cells_selected], s=2, norm=cnorm, cmap=plasma)
    # ax.set_xlim(umap_limits[0])
    # ax.set_ylim(umap_limits[1])
    # cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
    # cbar.set_label('Erreur de couverture')





    ################################################################################
    # list all the sub-rules with a null matching error (perfect match)
    # print('list of all the sub rules with perfect match:\n')
    # sub_df = df.loc[:,[instance.get_atom(atom_index)[0] for atom_index in LOLH_rule]]
    # sub_rules = []
    # for pos_sample in instance._pos_samples:
    #     sub_rule = []
    #     for atom_index in LOLH_rule:
    #         atom = instance.get_atom(atom_index)
    #         if sub_df[atom[0]][pos_sample] == atom[1]:
    #             sub_rule.append(atom_index)
    #     sub_rules.append(sub_rule)
    #     # print(pos_sample, ': ', [instance.get_atom(atom_index) for atom_index in sub_body])
    #
    # unique_sub, counts = np.unique(np.asarray(sub_rules), return_counts=True)
    # res = list(zip(counts, unique_sub))
    # res.sort(key=lambda elt:elt[0], reverse=True)
    # for elt in res:
    #     print(elt)
    #     print('\n')



    plt.show()

# NK_classification()
