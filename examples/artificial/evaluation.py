#!/usr/bin/python

import numpy as np
import pandas as pd
import random

# dimensionality reduction in the data
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib import cm


# LOLH objects
import sys
sys.path.append('../../python')

from instance import Instance
from solver import Solver
import histogram
import visualizer

# perform a comparision with PRIDE implementation of LFIT
import pylfit

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
def compute_A_B(instance, score):
    return (0, score*instance.n_negatives()), ((1.-score)*instance.n_positives(), instance.n_negatives())



#-------------------------------------------------------------------------------
def compute_rule_score(rule, instance):
    score = (0, 0)
    for atom_index in rule:
        atom_score = instance.atom_score[atom_index]
        score = (score[0] + atom_score[0]/instance.n_positives(), score[1] + atom_score[1]/instance.n_negatives())
    score = (score[0] / len(rule), score[1] / len(rule))
    return score


#-------------------------------------------------------------------------------
def dimensionality_reduction(df):

    X = df.values.copy()
    # standard scaling before PCA
    X=StandardScaler().fit_transform(X)
    # PCA with ten principal components
    X_pca = PCA(n_components=10).fit_transform(X)

    reducer = umap.UMAP(min_dist=0.3,n_neighbors=50,spread=1.0)
    embedding = reducer.fit_transform(X_pca)

    fig, axs = plt.subplots(1,2)

    col = ['forestgreen' for _ in range(int(df.shape[0]/2))] + ['darkred' for _ in range(int(df.shape[0]/2))]

    axs[0].scatter(embedding[:, 0], embedding[:, 1], s=0.5, c=col)
    # axs[0].set_aspect('equal', 'datalim')
    axs[0].set_title('projection UMAP')
    axs[0].set_xlabel('UMAP 1')
    axs[0].set_ylabel('UMAP 2')
    axs[0].legend(handles=[mpatches.Patch(color='forestgreen', label='exemples positifs'), mpatches.Patch(color='darkred', label='exemples négatifs')], loc='lower right')
    axs[0].set_aspect((axs[0].get_xlim()[1]-axs[0].get_xlim()[0])/(axs[0].get_ylim()[1]-axs[0].get_ylim()[0]))

    # PCA projection
    X = df.values.copy()
    X=StandardScaler().fit_transform(X)
    X_pca = PCA(n_components=2).fit_transform(X)
    axs[1].scatter(X_pca[:,0], X_pca[:,1], c=col, s=0.5)
    axs[1].set_title('projection PCA')
    axs[1].set_xlabel('PCA 1')
    axs[1].set_ylabel('PCA 2')
    axs[1].legend(handles=[mpatches.Patch(color='forestgreen', label='exemples positifs'), mpatches.Patch(color='darkred', label='exemples négatifs')])
    axs[1].set_aspect((axs[1].get_xlim()[1]-axs[1].get_xlim()[0])/(axs[1].get_ylim()[1]-axs[1].get_ylim()[0]))

    return

#-------------------------------------------------------------------------------
def plot_LOLH_PRIDE_scores(instance, LOLH_rule, LOLH_score, PRIDE_rules, PRIDE_scores):

    # plot the rule scores
    fig, ax = plt.subplots()

    ax.set_xlim((0, 1.))
    ax.set_ylim((0, 1.))
    ax.set_xlabel('score positif')
    ax.set_ylabel('score négatif')
    # ax.set_title('Comparaison des règles de PRIDE et LOLH')

    # diagonal indicating independant values
    ax.plot([0, 1], [0, 1], color='black', zorder=0)

    # plot the pride rules
    PRIDE_indexes = list(range(len(PRIDE_rules)))
    PRIDE_indexes.sort(key=lambda index: len(PRIDE_rules[index]))
    PRIDE_scores_sorted = [PRIDE_scores[index] for index in PRIDE_indexes] # sort all the points
    longest_PRIDE_rule = np.max([len(rule) for rule in PRIDE_rules])
    shortest_PRIDE_rule = np.min([len(rule) for rule in PRIDE_rules])
    colors = [len(PRIDE_rules[index]) for index in PRIDE_indexes]
    cnorm = mcolors.Normalize(vmin=shortest_PRIDE_rule, vmax=longest_PRIDE_rule)


    ax.scatter([score[0] for score in PRIDE_scores_sorted], [score[1] for score in PRIDE_scores_sorted], c=colors, norm=cnorm, cmap=plt.get_cmap('viridis'), marker='x',  zorder=1, label='règles PRIDE', facecolor='darkblue')
    cbar = fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap=plt.get_cmap('viridis')), ax=ax)
    cbar.set_label('longueur des règles PRIDE')
    ax.scatter(LOLH_score[0], LOLH_score[1], marker='*', s=70, color='chocolate', zorder=1, label='règle LOLH') # plot the optimized rule
    # plot the domination cone of the optimized rule
    ax.plot([LOLH_score[0], LOLH_score[0]], [LOLH_score[1], 0], linestyle='--', color='darkgrey', zorder=0)
    ax.plot([LOLH_score[0], 1.], [LOLH_score[1], LOLH_score[1]], linestyle='--', color='darkgrey', zorder=0)
    legend = ax.legend(loc='lower right')
    legend.legendHandles[0].set_color('darkcyan')

    return

#-------------------------------------------------------------------------------
def plot_atoms(instance, LOLH_atoms, PRIDE_atoms):

    # plot each atom score
    fig,ax = plt.subplots()

    ax.set_xlim((0,instance.n_positives()))
    ax.set_ylim((0,instance.n_negatives()))

    # ax.set_title('Atoms score')
    ax.set_xlabel('erreur positive')
    ax.set_ylabel('erreur négative')

    other_atoms = [atom_index for atom_index in range(instance.n_atoms()) if not atom_index in LOLH_atoms and not atom_index in PRIDE_atoms]

    cmap = plt.get_cmap('viridis')

    # plot LOLH atoms
    ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in LOLH_atoms], [instance.atom_score[ind_atom][1] for ind_atom in LOLH_atoms], marker='o', s=5, color='red', zorder=2, label='atomes LOLH')

    # plot pride atoms
    ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in PRIDE_atoms], [instance.atom_score[ind_atom][1] for ind_atom in PRIDE_atoms], marker='o', s=5, color='darkorange', zorder=2, label='atomes PRIDE')

    # other (not used) atoms
    ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in other_atoms], [instance.atom_score[ind_atom][1] for ind_atom in other_atoms], alpha=0.5, marker='x', zorder=0, label='autres atomes')

    # score = 0.0
    ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='black', label='score = 0.0', lw=1.5)

    # line for score = 0.6
    selection_score = 0.6
    # A = (0,selection_score*instance.n_negatives())
    # B = (instance.n_positives()*(1.-selection_score), instance.n_negatives())
    A,B = compute_A_B(instance, selection_score)
    ax.plot([A[0], B[0]], [A[1], B[1]], '--', color='black', label='score = 0.6', lw=1.5)

    ax.legend(loc='lower right')

    ax.set_aspect('equal')
    fig.tight_layout()

    return



#-------------------------------------------------------------------------------
def plot_LOLH_PRIDE_histograms(instance, LOLH_rule, PRIDE_rules, PRIDE_rule_scores):

    # histogram of the LOLH rule
    histo = histogram.Histogram(instance, LOLH_rule)
    fig, ax = plt.subplots()
    visualizer.plot_histograms(ax, histo, True)

    # selection of the best PRIDE rule
    PRIDE_indexes_sorted =  list(range(len(PRIDE_rules)))
    PRIDE_indexes_sorted.sort(key=lambda rule_index: PRIDE_rule_scores[rule_index][1] - PRIDE_rule_scores[rule_index][0], reverse=True)
    best_PRIDE_index = PRIDE_indexes_sorted[0]
    print('best pride rule: ', PRIDE_rules[best_PRIDE_index])

    # selection of the longest PRIDE rule
    longest_PRIDE_rule = PRIDE_rules[0]
    for PRIDE_rule in PRIDE_rules:
        if len(PRIDE_rule) > len(longest_PRIDE_rule):
            longest_PRIDE_rule = PRIDE_rule

    histo1 = histogram.Histogram(instance, longest_PRIDE_rule)
    histo2 = histogram.Histogram(instance, PRIDE_rules[best_PRIDE_index])

    # plot the histogram
    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1,0.9]})

    fig.tight_layout(w_pad=2)

    # ax.set_title('Histogrammes pour PRIDE')
    visualizer.plot_histograms(axs[0], histo1, True)
    axs[0].set_title(r'Plus longue règle PRIDE')

    visualizer.plot_histograms(axs[1], histo2, True)
    axs[1].set_title(r'Meilleure règle PRIDE')

    return


#-------------------------------------------------------------------------------
def multiobjective_comparison(instance, solver):

    body_length = 10
    biobj_score, biobj_bodies, biobj_weight = solver.compute_supported(body_length, 1, 1)
    biobj_score = [ (score[0]/float(instance.n_positives()), score[1]/float(instance.n_negatives())) for score in biobj_score]

    # selection of two "extreme" rules
    ind_left = 0
    for ind in range(len(biobj_score)):
        score = biobj_score[ind]
        if score[0] >= 0.4 and score[0] <= 0.6:
            if score[1] >= 6.4 and score[1] <= 6.8:
                print('left: ', ind)
                ind_left = ind
    print('left body: ', biobj_bodies[ind_left])
    left_score = biobj_score[ind_left]

    ind_right = 0
    for ind in range(len(biobj_score)):
        score = biobj_score[ind]
        if score[0] >= 3.6 and score[0] <= 4.0:
            if score[1] >= 9.5 and score[1] <= 9.7:
                print('right: ', ind)
                ind_right = ind
    print('right body: ', biobj_bodies[ind_right])
    right_score = biobj_score[ind_right]

    fig, ax = plt.subplots()
    # ax.set_title('Score des règles bi-objectives')
    ax.set_xlabel('score positif')
    ax.set_ylabel('score négatif')
    ind_reduced = [ind for ind in range(len(biobj_score)) if ind != ind_left and ind != ind_right]
    ax.scatter([biobj_score[ind][0] for ind in ind_reduced], [biobj_score[ind][1] for ind in ind_reduced], marker='x', s=25)
    ax.scatter([biobj_score[ind_left][0], biobj_score[ind_right][0]], [biobj_score[ind_left][1], biobj_score[ind_right][1]], color='red', marker='x', s=30)
    ax.text(5-0.3, 5+0.3, 'score négatif - score positif = 0', horizontalalignment='center', verticalalignment='center', rotation=45)
    ax.plot([0, body_length], [0, body_length], '--', color='k')

    ax.text(left_score[0] + 0.2, left_score[1]-0.2, r'$p_1$')
    ax.text(right_score[0] + 0.0, right_score[1]-0.6, r'$p_2$')

    delta = 0.5
    ax.set_xlim((0-delta, body_length+delta))
    ax.set_ylim((0-delta, body_length+delta))
    ax.set_aspect('equal')

    histo_left = histogram.Histogram(instance, biobj_bodies[ind_left])
    histo_right = histogram.Histogram(instance, biobj_bodies[ind_right])

    fig,axs = plt.subplots(1,2)
    fig.tight_layout(w_pad=3)
    axs[0].set_title(r'Histogramme de la règle $p_1$')
    visualizer.plot_histograms(axs[0], histo_left, True)
    axs[1].set_title(r'Histogramme de la règle $p_2$')
    visualizer.plot_histograms(axs[1], histo_right, True)

    return



#-------------------------------------------------------------------------------
def rule_histograms_comparisons(instance, solver):

    fig, ax = plt.subplots()
    ax.set_xlim([0, instance.n_positives()])
    ax.set_ylim([0, instance.n_negatives()])
    ax.scatter([instance.atom_score[ind][0] for ind in range(instance.n_atoms())], [instance.atom_score[ind][1] for ind in range(instance.n_atoms())], marker='x', zorder=1, alpha=0.5)
    # ax.set_title('Visualisation de plusieurs règles logiques')
    # ax.plot([0, nPos], [0, nNeg], '--', color='grey', zorder=3, alpha=0.7)
    ax.set_xlabel('erreur positive')
    ax.set_ylabel('erreur négative')

    fig, axsh = plt.subplots(2,2)
    fig.tight_layout(h_pad=4)


    score1 = 0.3
    score2 = 0.4
    A1, B1 = compute_A_B(instance, score1)
    A2, B2 = compute_A_B(instance, score2)
    selected_atoms, scores = solver.select_best_atoms_threshold(0.0)
    atoms_sandwich = [selected_atoms[ind] for ind in range(len(selected_atoms)) if scores[ind] >= score1 and scores[ind] <= score2]
    np.random.shuffle(atoms_sandwich)
    nAtomsSel = 10
    body_worst = atoms_sandwich[:nAtomsSel]
    print('body: ', body_worst)
    ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in body_worst], [instance.atom_score[ind_atom][1] for ind_atom in body_worst], marker='x', color='orangered', zorder=4)
    # histogramme
    histo = histogram.Histogram(instance, body_worst)
    axsh[0,0].set_title('Règle sous-optimale')
    visualizer.plot_histograms(axsh[0,0], histo, True)


    ax.fill([A1[0], B1[0], B2[0], A2[0]], [A1[1], B1[1], B2[1], A2[1]], zorder=3, color='green', alpha=0.3, label='règle sous-optimale')
    ax.plot([A1[0], B1[0]], [A1[1], B1[1]], '--', color='k', zorder=3)
    ax.plot([A2[0], B2[0]], [A2[1], B2[1]], '--', color='k', zorder=3)

    score1 = -0.7
    score2 = -0.6
    A1, B1 = compute_A_B(instance, score1)

    A2, B2 = compute_A_B(instance, score2)

    selected_atoms, scores = solver.select_best_atoms_threshold(-1.)
    atoms_sandwich = [selected_atoms[ind] for ind in range(len(selected_atoms)) if scores[ind] >= score1 and scores[ind] <= score2]
    np.random.shuffle(atoms_sandwich)
    nAtomsSel = 10
    body_worst = atoms_sandwich[:nAtomsSel]
    print('body: ', body_worst)
    ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in body_worst], [instance.atom_score[ind_atom][1] for ind_atom in body_worst], marker='x', color='orangered', zorder=4)
    # histogramme
    histo = histogram.Histogram(instance, body_worst)
    axsh[0,1].set_title('Règle inverse')
    visualizer.plot_histograms(axsh[0,1], histo, True)

    ax.fill([A1[0], B1[0], B2[0], A2[0]], [A1[1], B1[1], B2[1], A2[1]], zorder=3, color='orange', alpha=0.3, label='règle inverse')
    ax.plot([A1[0], B1[0]], [A1[1], B1[1]], '--', color='k', zorder=3)
    ax.plot([A2[0], B2[0]], [A2[1], B2[1]], '--', color='k', zorder=3)

    score1 = -0.05
    score2 = 0.05
    A1, B1 = compute_A_B(instance, score1)
    A2, B2 = compute_A_B(instance, score2)
    selected_atoms, scores = solver.select_best_atoms_threshold(-1.)
    atoms_sandwich = [selected_atoms[ind] for ind in range(len(selected_atoms)) if scores[ind] >= score1 and scores[ind] <= score2]
    np.random.shuffle(atoms_sandwich)
    nAtomsSel = 10
    body_worst = atoms_sandwich[:nAtomsSel]
    print('body: ', body_worst)
    ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in body_worst], [instance.atom_score[ind_atom][1] for ind_atom in body_worst], marker='x', color='orangered', zorder=4, label='atomes sélectionnés')
    # histogramme
    histo = histogram.Histogram(instance, body_worst)
    axsh[1,0].set_title('Règle inappropriée')
    visualizer.plot_histograms(axsh[1,0], histo, True)


    ax.fill([A1[0], B1[0], B2[0], A2[0]], [A1[1], B1[1], B2[1], A2[1]], zorder=3, color='red', alpha=0.3, label='règle inappropriée')
    ax.plot([A1[0], B1[0]], [A1[1], B1[1]], '--', color='k', zorder=3)
    ax.plot([A2[0], B2[0]], [A2[1], B2[1]], '--', color='k', zorder=3)

    atomes = [ind for ind in range(instance.n_atoms())]
    np.random.shuffle(atomes)
    body_random = atomes[:10]
    histo = histogram.Histogram(instance, body_random)
    axsh[1,1].set_title('Règle aléatoire')
    visualizer.plot_histograms(axsh[1,1], histo, True)


    legend = ax.legend(loc='lower left')

    handles_old, labels_old = ax.get_legend_handles_labels()

    new_handles = []
    for ind_handle in range(len(handles_old)-1):
        patch = handles_old[ind_handle]
        new_handles.append(Patch(edgecolor='k', facecolor=colorConverter.to_rgba(patch.get_facecolor(), alpha=0.3), linewidth=0.5, label=labels_old[ind_handle]))
    new_handles.append(handles_old[-1])
    ax.legend(loc='lower left', handles=new_handles, labels=labels_old)

    return

#-------------------------------------------------------------------------------
def compute_matching_error(df, instance, state, rule):
    matching_error = 0
    for atom_index in rule:
        atom = instance.get_atom(atom_index)
        if df[atom[0]][state] != atom[1]:
            matching_error += 1
    return matching_error

#-------------------------------------------------------------------------------
def k_fold_cross_validation(df, instance, k):

    # create nsamples/k subset of the data
    pos_indexes = instance._pos_samples
    pos_sub_indexes = [[] for _ in range(k)]
    for ind in range(len(pos_indexes)):
        k_prime = int((ind/instance.n_positives())*k)
        pos_sub_indexes[k_prime].append(pos_indexes[ind])

    neg_indexes = instance._neg_samples
    neg_sub_indexes = [[] for _ in range(k)]
    for ind in range(len(neg_indexes)):
        k_prime = int((ind/instance.n_negatives())*k)
        neg_sub_indexes[k_prime].append(neg_indexes[ind])

    LOLH_perf = [[],[]]
    PRIDE_perf = [[],[]]
    RANDOM_perf = [[],[]]

    # for ind_k in range(1):

    for ind_k in range(k):

        # create the training and test set
        #ind_K = 0
        train_pos = pos_sub_indexes[ind_k]
        train_neg = neg_sub_indexes[ind_k]
        test_pos = []
        test_neg = []
        for ind in range(k):
            if ind != ind_k:
                test_pos += pos_sub_indexes[ind]
                test_neg += neg_sub_indexes[ind]

        # debug
        # test_pos = train_pos
        # test_neg = train_neg

        df_train = df.loc[train_pos+train_neg]
        df_test = df.loc[test_pos+test_neg]

        # training set
        inst_train = Instance.create_instance_explicit(df_train, train_pos, train_neg)

         # test set
        inst_test = Instance.create_instance_explicit(df_test, test_pos, test_neg)

        # compute the model on the training set
        solver = Solver(inst_train)

        # computation of the solution with LOLH and PRIDE
        # the atom indexes are converted from the train instance to the test instance

        # compute LOLH body
        # LOLH_body, scores = solver.select_best_atoms_threshold(0.75)
        LOLH_body, scores = solver.select_best_atoms_threshold(0.6)
        LOLH_body = [inst_train.get_atom(atom_index) for atom_index in LOLH_body]
        LOLH_body = [inst_test.get_atom_index(atom) for atom in LOLH_body]

        # compute PRIDE bodies
        PRIDE_bodies = compute_pride_bodies(inst_train)
        PRIDE_bodies = [[inst_train.get_atom(atom_index) for atom_index in pride_body] for pride_body in PRIDE_bodies]
        PRIDE_bodies = [[inst_test.get_atom_index(atom) for atom in PRIDE_body] for PRIDE_body in PRIDE_bodies]

        # compute a random rule
        RANDOM_rule = list(range(inst_train.n_atoms()))
        np.random.shuffle(RANDOM_rule)
        RANDOM_rule = RANDOM_rule[:len(LOLH_body)]
        # print('random rule: ', RANDOM_rule)

        # histogram of the random rule on the test instance
        # histo = histogram.Histogram(inst_test, RANDOM_rule)
        # fig, ax = plt.subplots()
        # ax.set_title('random histogram k=' + str(ind_k))
        # visualizer.plot_histograms(ax, histo, True)

        # error of the atoms
        # fig, ax = plt.subplots()
        # ax.scatter([inst_test.atom_score[ind_atom][0] for ind_atom in range(inst_test.n_atoms())], [inst_test.atom_score[ind_atom][1] for ind_atom in range(inst_test.n_atoms())], marker='x')
        # ax.scatter([inst_test.atom_score[ind_atom][0] for ind_atom in RANDOM_rule], [inst_test.atom_score[ind_atom][1] for ind_atom in RANDOM_rule], marker='x', c='red')
        # ax.plot([0, inst_test.n_positives()], [0, inst_test.n_negatives()], color='k')

        # atoms in the train and test instance
        atoms_train_test = [inst_train.get_atom(atom_index) for atom_index in range(inst_train.n_atoms()) if inst_test.has_atom(inst_train.get_atom(atom_index))]
        col = []
        for atom in atoms_train_test:
            ind_atom_train = inst_train.get_atom_index(atom)
            score = -inst_train.atom_score[ind_atom_train][0]*len(train_neg)+inst_train.atom_score[ind_atom_train][1]*len(train_pos)
            col.append(score)

        # LOLH performance
        pos_match = 0
        neg_match = 0
        for pos in test_pos:
            matching_error = compute_matching_error(df_test, inst_test, pos, LOLH_body)
            if matching_error < len(LOLH_body)/2.:
                pos_match += 1
        for neg in test_neg:
            matching_error = compute_matching_error(df_test, inst_test, neg, LOLH_body)
            if matching_error < len(LOLH_body)/2.:
                neg_match += 1

        # print('pos match: ', pos_match)
        # print('neg match: ', neg_match)
        LOLH_perf[0].append(float(pos_match)/len(test_pos))
        LOLH_perf[1].append(float(neg_match)/len(test_neg))

        # PRIDE performance
        pos_match = 0
        neg_match = 0

        for pos in test_pos:
            for pride_body in PRIDE_bodies:
                matching_error = compute_matching_error(df_test, inst_test, pos, pride_body)
                if matching_error == 0:
                    pos_match += 1
                    break
        for neg in test_neg:
            for pride_body in PRIDE_bodies:
                matching_error = compute_matching_error(df_test, inst_test, neg, pride_body)
                if matching_error == 0:
                    neg_match += 1
                    break

        PRIDE_perf[0].append(pos_match/len(test_pos))
        PRIDE_perf[1].append(neg_match/len(test_neg))

        # RANDOM model performance
        pos_match = 0
        neg_match = 0
        for pos in test_pos:
            matching_error = compute_matching_error(df_test, inst_test, pos, RANDOM_rule)
            if matching_error < len(RANDOM_rule)/2.:
                pos_match += 1
        for neg in test_neg:
            matching_error = compute_matching_error(df_test, inst_test, neg, RANDOM_rule)
            if matching_error < len(RANDOM_rule)/2.:
                neg_match += 1

        # print('pos match: ', pos_match)
        # print('neg match: ', neg_match)

        RANDOM_perf[0].append(float(pos_match)/len(test_pos))
        RANDOM_perf[1].append(float(neg_match)/len(test_neg))

    print('Performances of LOLH:\n')
    print(LOLH_perf)
    print('\n')

    print('Performances of PRIDE:\n')
    print(PRIDE_perf)
    print('\n')

    print('Performances of a RANDOM rule:\n')
    print(RANDOM_perf)

    return

#-------------------------------------------------------------------------------
def evaluate(filename):

    random_seed = 42
    random.seed(random_seed)
    np.random.seed(random_seed)

    # load the matrix
    df = pd.read_csv(filename, index_col=0)




    ########################################################
    # dimensionality reduction to visualize the data
    print('*------------------------------------------------------------*')
    print('    Dimensionality reduction on the data for visualization    ')
    print('*------------------------------------------------------------*')

    dimensionality_reduction(df)

    plt.show()

    print('\n\n\n\n')




    # classification from positive vs negative
    pos_samples = ['s_'+str(ind) for ind in range(0, int(df.shape[0]/2.))]
    neg_samples = ['s_'+str(ind) for ind in range(int(df.shape[0]/2.)+1, df.shape[0])]

    instance_name = 'artificial'

    inst = Instance.create_instance_explicit(df, pos_samples, neg_samples)

    # extract positive/negative solutions for solution visualisation
    positive_cells = inst._pos_samples
    negative_cells = inst._neg_samples

    print('Original instance:')
    # print('n atoms: ', inst.n_atoms())
    print('n positive samples: ', inst.n_positives())
    print('n negative samples: ', inst.n_negatives())

    # LOLH solver

    solver = Solver(inst)



    ########################################################
    # computation of LOLH rule

    # selection score = 0.75
    selection_score = 0.6

    selected_atoms, scores = solver.select_best_atoms_threshold(selection_score)

    print('n atoms selected: ', len(selected_atoms))

    LOLH_rule = selected_atoms

    # compute the rule score for the LOLH rule
    LOLH_rule_score = compute_rule_score(LOLH_rule, inst)

    # plot the body
    # print([inst.get_atom(ind) for ind in LOLH_rule])




    ########################################################
    # computation of PRIDE rules
    PRIDE_rules = compute_pride_bodies(inst)
    PRIDE_atoms = list(np.unique([atom_index for pride_rule in PRIDE_rules for atom_index in pride_rule]))
    print('n rules pride: ', len(PRIDE_rules))
    print('rule 0: ', PRIDE_rules[0])

    # compute a global score for each rule
    PRIDE_rule_scores = []
    for PRIDE_rule in PRIDE_rules:
        PRIDE_rule_scores.append(compute_rule_score(PRIDE_rule, inst))





    ########################################################
    # plot LOLH and PRIDE rule scores

    print('*-------------------------------------------------------*')
    print('    Comparison of the rule scores from LOLH and PRIDE    ')
    print('*-------------------------------------------------------*')

    plot_LOLH_PRIDE_scores(inst, LOLH_rule, LOLH_rule_score, PRIDE_rules, PRIDE_rule_scores)

    plt.show()

    print('\n\n\n\n')



    ########################################################
    # plot all the atoms
    print('*-------------------------------------------------*')
    print('    Comparison of the atoms from LOLH and PRIDE    ')
    print('*-------------------------------------------------*')

    plot_atoms(inst, LOLH_rule, PRIDE_atoms)

    plt.show()

    print('\n\n\n\n')




    ########################################################
    # plot the histograms of the LOLH rule and of some PRIDE rules
    print('*------------------------------------------------------*')
    print('    Comparison of the histograms from LOLH and PRIDE    ')
    print('*------------------------------------------------------*')

    plot_LOLH_PRIDE_histograms(inst, LOLH_rule, PRIDE_rules, PRIDE_rule_scores)

    plt.show()

    print('\n\n\n\n')



    ################################################################
    # computation of the multi-objective rules
    print('*-------------------------------------------------------*')
    print('    Comparison of two different multi-objective rules    ')
    print('*-------------------------------------------------------*')

    multiobjective_comparison(inst, solver)

    plt.show()

    print('\n\n\n\n')




    ################################################################
    # display several rules with histograms
    print('*-----------------------------------------------------*')
    print('    Comparison of the histograms from several rules    ')
    print('*-----------------------------------------------------*')

    rule_histograms_comparisons(inst, solver)

    plt.show()

    print('\n\n\n\n')




    ################################################################
    # K-fold validation on PRIDE and LOLH
    print('*----------------------------------------------------------------*')
    print('    Comparison of several models with 10-fold cross validation    ')
    print('*----------------------------------------------------------------*')


    k = 10
    k_fold_cross_validation(df, inst, k)

    plt.show()

    print('\n\n\n\n')

    return

# evaluate('../../dataset/artificial/artificial_matrix.csv')
