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

import pylfit # perform a comparision with PRIDE implementation of LFIT


#------------------------------------------------------------------------------
def produce_visuals(normalized_dataset, instance, body):

    # create an histogram
    histo = histogram.Histogram(instance, body)
    fig, ax = plt.subplots()
    visualizer.plot_histograms(ax, histo)
    ax.set_title('LOLH rule matching error for the NK cells')
    ax.set_ylim([0, 0.3])

    # plot the violins
    fig, ax = plt.subplots()
    visualizer.plot_violins(ax, normalized_dataset, instance, body)
    ax.set_title('Violins of the selected genes from LOLH (NK cells)')

    return


#------------------------------------------------------------------------------
def export_genes(instance, sorted_atoms, scores, threshold):

    # export the selected genes in a file

    # export the genes in a file
    file = open('../../dataset/Imagine/cell_classification/NK_gene_selection.csv', 'w')
    file.write('gene,discrete_value,score\n')
    ind = 0
    while scores[ind] >= threshold:
        atom = instance.get_atom(sorted_atoms[ind])
        file.write(atom[0] + ',' + str(atom[1]) + ',' + str(scores[ind])+'\n')
        ind += 1
    file.close()

    return

#------------------------------------------------------------------------------
def compute_pride_bodies(instance):

    # computation of the bodies of all rules returned from pride algorithm on the instance

    variable = -1
    value = -1 # not important, just for rules'head

    nb_features = instance.dataset.shape[0] # nb columns in the dataframe

    positives = instance.dataset.loc[instance._pos_samples].values.tolist()
    negatives = instance.dataset.loc[instance._neg_samples].values.tolist()

    verbose = 0
    output = pylfit.algorithms.pride.PRIDE.fit_var_val(variable, value, nb_features, positives, negatives, verbose)

    # for rule in output:
    #    print(rule)
    # print( str(len(output)) + ' rules induced by PRIDE')



    # print(rule._body_variables)
    # to_disp = ''
    # for var_index in rule._body_variables:
    #     to_disp += str(var_index) + '_' + str(rule._body_values[var_index]) + ', '
    # print(to_disp)

    # compute the bodies of PRIDE rules: indexes of the selected atoms from the instance of NK cells
    pride_bodies = []
    for rule in output:

        pride_body = []

        for var_index in rule._body_variables:
            gene = instance.dataset.columns[var_index]
            value = rule._body_values[var_index]
            atom_index = instance.get_atom_index((gene, value))
            pride_body.append(atom_index)

        pride_bodies.append(pride_body)

    return pride_bodies

#------------------------------------------------------------------------------
def compare_pride(instance, lolh_body, pride_bodies):

    # compare the genes selected by LOLH with the rule computed with pride

    # compute a global score for each rule
    pride_rule_scores = []

    for pride_body in pride_bodies:

        score = (0, 0)

        for atom_index in pride_body:
            atom_score = instance.atom_score[atom_index]
            score = (score[0] + atom_score[0], score[1] + atom_score[1])

        score = (score[0] / len(pride_body), score[1] / len(pride_body))
        pride_rule_scores.append(score)


    # compute the rule score for the LOLH rule
    lolh_rule_score = (0,0)
    for index in lolh_body:
        atom_score = instance.atom_score[index]
        lolh_rule_score = (lolh_rule_score[0] + atom_score[0], lolh_rule_score[1] + atom_score[1])
    lolh_rule_score = (lolh_rule_score[0]/len(lolh_body), lolh_rule_score[1]/len(lolh_body))


    # display the normalized scores for all rules (between (0,0) and (n_positives, n_negatives))
    fig, ax = plt.subplots()

    ax.set_xlim((0, instance.n_positives()))
    ax.set_ylim((0, instance.n_negatives()))
    ax.set_xlabel('rule positive score')
    ax.set_ylabel('rule negative score')
    ax.set_title('Comparison of PRIDE and LOLH rule scores')


    ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='black', zorder=0) # diagonal indicating independant values

    # plot the pride rules
    pride_indexes = [ind for ind in range(len(pride_bodies))]
    pride_indexes.sort(key=lambda ind: len(pride_bodies[ind]))
    pride_rule_scores_sorted = [pride_rule_scores[ind] for ind in pride_indexes] # sort all the points
    longest_rule = np.max([len(body) for body in pride_bodies])
    colors = [len(pride_bodies[ind]) for ind in pride_indexes]
    cnorm = mcolors.Normalize(vmin=0, vmax=longest_rule)

    ax.scatter([score[0] for score in pride_rule_scores_sorted], [score[1] for score in pride_rule_scores_sorted], c=colors, norm=cnorm, cmap=plt.get_cmap('viridis'), marker='x',  zorder=1, label='PRIDE rules', facecolor='darkblue')
    cbar = fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap=plt.get_cmap('viridis')), ax=ax)
    cbar.set_label('PRIDE rule body length')
    ax.scatter(lolh_rule_score[0], lolh_rule_score[1], marker='*', s=70, color='chocolate', zorder=1, label='LOLH rule') # plot the optimized rule
    # plot the domination cone of the optimized rule
    ax.plot([lolh_rule_score[0], lolh_rule_score[0]], [lolh_rule_score[1], 0], linestyle='--', color='darkgrey', zorder=0)
    ax.plot([lolh_rule_score[0], instance.n_positives()], [lolh_rule_score[1], lolh_rule_score[1]], linestyle='--', color='darkgrey', zorder=0)
    legend = ax.legend(loc='lower right')
    legend.legendHandles[0].set_color('darkcyan')


    # plot the individual atom scores, with color to distinguish from PRIDE and optimized rule
    pride_atoms = list(np.unique(np.concatenate(pride_bodies)))
    fig, ax = plt.subplots()
    ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='k', zorder=1) # baseline
    ax.scatter([instance.atom_score[index][0] for index in lolh_body], [instance.atom_score[index][1] for index in lolh_body], marker='o', facecolor='darkorange', edgecolor='black', linewidth=0.5, zorder=2, label='LOLH atoms')
    ax.scatter([instance.atom_score[index][0] for index in pride_atoms], [instance.atom_score[index][1] for index in pride_atoms], marker='s', facecolor='dodgerblue', edgecolor='black', linewidth=0.5, zorder=2, label='PRIDE atoms')
    ax.scatter([instance.atom_score[index][0] for index in range(instance.n_atoms()) if not index in pride_atoms+lolh_body], [instance.atom_score[index][1] for index in range(instance.n_atoms()) if not index in pride_atoms+lolh_body], color='grey', marker='x', zorder=0, label='other atoms')
    ax.legend(loc='lower right')
    ax.set_xlabel('atom positive error')
    ax.set_ylabel('atom negative error')
    ax.set_title('Individual error of the atoms from PRIDE and LOLH rules')

    return


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

    # read the normalized matrix
    filename = '../../dataset/Imagine/normalized_matrix.csv'
    df_normalized = pd.read_csv(filename, index_col=0).T
    # print(df_normalized.head())

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
    # print('threshold: ', threshold)
    # print(scores[:10])
    threshold *= instance.n_positives()*instance.n_negatives()
    lolh_body = [sorted_atoms[ind] for ind in range(len(sorted_atoms)) if scores[ind] >= threshold] # 21 atoms selected

    print('-> Output of LOLH: ', [instance.get_atom(index)[0] + '_' + str(instance.get_atom(index)[1]) for index in lolh_body])

    print('Exportation of the genes selected by PRIDE \n')

    # compute the bodies from pride rules
    pride_bodies = compute_pride_bodies(instance)

    # compute the different atoms in the PRIDE rules
    pride_atoms = list(np.unique([instance.get_atom(atom_index) for body in pride_bodies for atom_index in body]))

    # export the genes selected by PRIDE in a txt file
    pride_genes = list(np.unique([instance.get_atom(atom_index)[0] for body in pride_bodies for atom_index in body]))
    file = open('../../dataset/Imagine/cell_classification/PRIDE_genes.txt', 'w')
    for index in range(len(pride_genes)):
        file.write(pride_genes[index])
        if index < len(pride_genes)-1:
            file.write(',')
    file.write('\n')
    file.close()

    print('\n')

    print('-> PRIDE algorithm returned ', len(pride_bodies), ' rules, using in total ', len(pride_atoms), ' atoms, and ', len(pride_genes), ' genes')

    print('\n')


    print('Visualisations of the LOLH rule')

    print('\n')

    # Visualisation of the LOLH classification rule
    produce_visuals(df_normalized, instance, lolh_body)

    print('Comparison between PRIDE and LOLH')

    print('\n')

    # perform comparisions between LOLH and PRIDE
    compare_pride(instance, lolh_body, pride_bodies)

    plt.show()

# celltype_classification()
