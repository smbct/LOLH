#!/usr/bin/python3

import numpy as np
import pandas as pd

import random

from instance import Instance
from solver import Solver

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches

from matplotlib import cm
from matplotlib import colors as mcolors


import visualizer
import histogram

import math

from numpy.random import default_rng

from scipy.stats import t
from scipy.stats import levy_stable
from scipy.stats import cauchy

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap

import pylfit # perform a comparision with PRIDE implementation of LFIT

#---------------------------------------------------------------------
def singlecell_matrix_experiment():
    
    random.seed(42)
    np.random.seed(42)
    
    ####################################################################
    ###### read cell types
    ####################################################################
    file_name = '../dataset/IMAGINE/cell_types.csv'
    df_cell_types = pd.read_csv(file_name, index_col=0)
    df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)
    print(df_cell_types.head())
    
    ####################################################################
    ###### read macro cell types
    ####################################################################
    file_name = '../dataset/IMAGINE/cell_types_macro.csv'
    df_cell_types_macro = pd.read_csv(file_name, index_col=0)
    df_cell_types_macro.rename(columns={'cellType_macro': 'Label'}, inplace=True)
    print(df_cell_types_macro.head())
    
    
    ####################################################################
    ###### read the normalized dataset
    ####################################################################
    file_name = '../../../Learning/IMAGINE_dataset/dataset/ctls_normalised_counts.csv'
    
    df_normalised = pd.read_csv(file_name, index_col = 0).T
    print(df_normalised.head())
    
    ####################################################################
    ###### read the discretized dataset
    ####################################################################
    
    file_name = '../dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv'
    
    # Load the discrete dataset
    df = pd.read_csv(file_name, index_col=0)
        
    # select only the T cells
    #df = df.loc[df_cell_types_macro.loc[df_cell_types_macro['Label']=='T'].index]
    
    # remove axis name (name of the previously first column)
    df.rename_axis('Barcode', axis=0, inplace=True)
    
    # select a sub dataset
    # df = df.loc[cd4_negatives+cd4_positives]
    
    print(df.head())
    
    
    
    # classification of CD8+
    inst = Instance.create_cluster_instance(df, df_cell_types, 'NK')
    #inst = Instance.create_cluster_instance(df, df_cell_types, 'CD14')
    instance_name = 'NK'
    
    # inst = Instance.create_random_instance(df, 0.3)
    
    # extract positive/negative solutions for solution visualisation
    positive_cells = inst._pos_samples
    negative_cells = inst._neg_samples
    
    print('Original instance:')
    # print('n atoms: ', inst.n_atoms())
    print('n positive samples: ', inst.n_positives())
    print('n negative samples: ', inst.n_negatives())
    
    solver = Solver(inst)
    
    
    ############################################
    ### computation of supported front
    ############################################
    
    body_length = 20
    
    # create line equations from these points
    x1, y1 = 0.,0.
    x2, y2 = (inst.n_positives(), inst.n_negatives())
    a = 1.
    b = (x1-x2)/(y2-y1)
    c = -a*x1 - b*y1
    
    max_dist = a*0+b*y2+c
    
    # compute a line at distance 0.3
    p031 = [0,0]
    p031[1] = (0.3*max_dist-c)/b
    
    p032 = [0,y2]
    p032[0] = (0.3*max_dist-c-b*p032[1])/a
    
    # compute the list of distances to the diagonal
    atom_distance = [-1 for ind_atom in range(inst.n_atoms())]
    for ind_atom in range(inst.n_atoms()):
        score = inst.atom_score[ind_atom]
        atom_distance[ind_atom] = (a*score[0] + b*score[1] + c)/max_dist
    sorted_atoms = [ind for ind in range(inst.n_atoms())]
    sorted_atoms.sort(key=lambda ind:atom_distance[ind], reverse=True)
    
    
    stop_ind = 0
    for elt in sorted_atoms:
        if atom_distance[elt] < 0.3:
            break
        else:
            stop_ind += 1
    
    sorted_atoms = sorted_atoms[:stop_ind]
    
    #selected_atoms, scores = solver.select_k_best_atoms(50)
    selected_atoms, scores = solver.select_best_atoms_threshold(0.3)
    
    # print([(inst.get_atom(ind), atom_distance[ind]) for ind in sorted_atoms])
    
    print('n atoms selected: ', len(selected_atoms))
    #for elt in [(inst.get_atom(ind), atom_distance[ind]) for ind in selected_atoms]:
    #    print(elt)
    
    body_length = 20
    body = selected_atoms[:body_length]
    
    # plot the body
    print([inst.get_atom(ind) for ind in body])
    
    # create an histogram to visualize rule error over samples
    histo = histogram.Histogram(inst, body)
    
    
    # plot the histogram
    fig, ax = plt.subplots()
    ax.set_title('Histograms for ' + instance_name + ' classification')
    visualizer.plot_histograms(ax, histo, True)
    
    nPos = inst.n_positives()
    nNeg = inst.n_negatives()
    
    
    # plot each atom score
    fig,ax = plt.subplots()
    
    ax.set_xlim((0,nPos))
    ax.set_ylim((0,nNeg))
    
    #ax.set_aspect('equal')
    
    # ax.set_title('Atoms score')
    # ax.set_xlabel('Positive score')
    # ax.set_ylabel('Negative score')
    
    ax.set_title('Erreurs des atomes logiques (classification de NK)')
    ax.set_xlabel('Erreur positive')
    ax.set_ylabel('Erreur négative')
    
    
    
    # colors = ['royalblue' for _ in range(len(inst.atom_score))]
    ax.scatter([elt[0] for elt in inst.atom_score], [elt[1] for elt in inst.atom_score], marker='x')
    
    # y / x
    # ax.plot([0, inst.n_positives()], [0, inst.n_negatives()], color='red')
    
    # 0.3 threshold
    # ax.plot([p031[0], p032[0]], [p031[1], p032[1]], color='lightseagreen')
    
    # compute atom distance to the diagonal vs its lateral position on the diagonal
    nPos = inst.n_positives()
    nNeg = inst.n_negatives()
    
    a_diag = nNeg
    b_diag = -nPos
    c_diag = 0
    
    
    # diagonal
    ax.plot([0, nPos], [0, -a_diag*nPos/b_diag], color='red')
    
    ax.legend(['score = 0.0'], loc='lower right')
    
    # orthogonal projection
    x_mid = nPos/2
    y_mid = nNeg/2
    
    # compute the relative distances do diagonal and normal for each atoms
    atom_relative_distances = [(0,0) for _ in range(inst.n_atoms())]
    
    max_dist = -a_diag*0 -b_diag*nNeg
    
    k = math.sqrt(a_diag**2 + b_diag**2)
    for atom_index in range(inst.n_atoms()):
        pos = inst.atom_score[atom_index]
        c = -a_diag*pos[0] - b_diag*pos[1]
        dist_diag = c/max_dist
        if c > 0:
            x1 = 0
            y1 = (-a_diag*x1-c)/b_diag
            y2 = nNeg
            x2 = (-b_diag*y2-c)/a_diag
        else:
            x2 = nPos
            y2 = (-a_diag*x2-c)/b_diag
            y1 = 0
            x1 = (-b_diag*y1-c)/a_diag
            
        dist_left = math.sqrt((pos[0]-x1)**2 + (pos[1]-y1)**2)
        length = math.sqrt((x2-x1)**2 + (y2-y1)**2)
        lateral_dist = dist_left/length
    
        atom_relative_distances[atom_index] = (dist_diag, lateral_dist)
    
    colors = [max(inst.atom_score[ind_atom][0],inst.atom_score[ind_atom][1])/(inst.atom_score[ind_atom][0]+inst.atom_score[ind_atom][1]) for ind_atom in range(inst.n_atoms())]
    
    fig, ax = plt.subplots()
    ax.scatter([elt[0] for elt in atom_relative_distances], [elt[1] for elt in atom_relative_distances], marker='x', c=colors)
    ax.set_xlabel('distance to diagonal')
    ax.set_ylabel('lateral position')
    
    # add histograms
    
    # verif: atom with colors
    # plot each atom score
    fig,ax = plt.subplots()
    # ax.set_aspect('equal')
    ax.set_xlim((0,nPos))
    ax.set_ylim((0,nNeg))
    ax.set_title('Atoms score')
    ax.set_xlabel('Positive score')
    ax.set_ylabel('Negative score')

    
    ax.scatter([elt[0] for elt in inst.atom_score], [elt[1] for elt in inst.atom_score], marker='x', c=[elt[1] for elt in atom_relative_distances])
    
    ax.plot([0, nPos], [0, nNeg], color='red')
    
    selected_atoms = [ind for ind in range(inst.n_atoms())]
    selected_atoms.sort(key=lambda ind: atom_relative_distances[ind][0], reverse=False)
    
    fig, axs = plt.subplots(3)
    axs[0].hist([elt[0] for elt in atom_relative_distances], edgecolor='k', bins=100)
    axs[1].hist([elt[1] for elt in atom_relative_distances], edgecolor='k', bins=100)
    axs[2].hist([max(inst.atom_score[ind_atom][0],inst.atom_score[ind_atom][1])/(inst.atom_score[ind_atom][0]+inst.atom_score[ind_atom][1]) for atom_index in range(inst.n_atoms())], edgecolor='k', bins=100)
    
    plt.show()




###############################################################
def artificial_data_experiment():
    
    nPos = 1000
    nNeg = 1000
    
    nVar = 10000
    nAtoms = nVar*2
    
    # random_seed = 46
    random_seed = 52 # worst dataset
    
    random.seed(random_seed)
    np.random.seed(random_seed)
    
    # sampling over the normal distribution
    
    
    rng = default_rng(random_seed)
    
    # mu, sigma = 0, 0.1 # mean and standard deviation
    # sample_normal = rng.normal(mu, sigma, nAtoms)
    
    sample_uniform = rng.uniform(0, 1, nVar*2)
    
    df = 10
    sample_score = t.rvs(df, size=nVar)
    # sample_student = [elt for elt in sample_student if elt >= 0]

    # https://en.wikipedia.org/wiki/Stable_distribution
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.levy_stable.html
    # stable distribution
    # alpha, beta = 1.995, 0.
    # sample_score = levy_stable.rvs(alpha, beta, size=nVar)
    
    # x0 = 0
    # gamma = 0.5
    # global sample_score
    # sample_score = cauchy.rvs(loc=x0, scale=gamma, size=nVar)
    
    # normal is normalized
    max_normal = np.max(sample_score)
    min_normal = np.min(sample_score)
    
    print(min_normal, max_normal)
    
    max_abs = max(max_normal, -min_normal) + 0.2
    for ind in range(len(sample_score)):
        sample_score[ind] = sample_score[ind]/max_abs
        sample_score[ind] *= 0.70
    # threshold = 50
    # for ind in range(len(sample_score)):
    #     sample_score[ind] = sample_score[ind]/threshold
    
    fig, axs = plt.subplots(2)
    axs[0].hist(sample_score, edgecolor='k', bins=100)
    axs[1].hist(sample_uniform, edgecolor='k', bins=100)

    
    #print(sample_normal)

    # fig,ax = plt.subplots()
    # ax.scatter(sample_normal, sample_uniform, marker='x')
    # ax.set_xlim([-1, 1])
    
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
        
        
    # generating positiv and negative examples
    global positives
    global negatives
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
        
        
    atom_scores_sym = []
    for score in atom_scores:
        atom_scores_sym.append( (nPos-score[0], nNeg-score[1]) )
    atom_scores += atom_scores_sym

    fig,ax = plt.subplots()
    ax.scatter([elt[0] for elt in atom_scores], [elt[1] for elt in atom_scores], marker='x')
    ax.plot([0, nPos], [0, nNeg], color='red')    
    ax.set_xlim([0, nPos])
    ax.set_ylim([0, nNeg])
    
    ind_atoms_sorted = [ind for ind in range(len(atom_scores))]
    score = [-nNeg*error[0] +nPos*error[1] for error in atom_scores]
    ind_atoms_sorted.sort(key = lambda ind: score[ind], reverse=True)
    fig, ax = plt.subplots()
    ax.scatter([atom_scores[ind][0] for ind in ind_atoms_sorted], [atom_scores[ind][1] for ind in ind_atoms_sorted], c=[score[ind] for ind in ind_atoms_sorted], marker='x')
    ax.set_xlim([0, nPos])
    ax.set_ylim([0, nNeg])
    ax.plot([0, nPos], [0, nNeg], color='red')
    
    selected_atoms = ind_atoms_sorted[:20]
    
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
        pos_histogram[score] /= float(nPos)
        
    for score in range(len(neg_histogram)):
        neg_histogram[score] /= float(nNeg)
    
    print('positive histogram: ')
    print(pos_histogram)
    
    print('negative histogram:')
    print(neg_histogram)
    
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
    
    
    # create a csv file with the matrix (ordered positive and negative examples)
    matrix = positives + negatives
    variables = ['v_'+str(ind) for ind in range(nVar)]
    examples = ['s_' + str(ind) for ind in range(nPos+nNeg)]
    
    # exportation into a dataframe
    # global dataframe
    # dataframe = pd.DataFrame(matrix, columns = variables, index= examples)
    # dataframe.to_csv('../dataset/artificial_matrix_worst.csv')
    
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
def LOLH_artificial_data():
    
    
    
    random.seed(42)
    np.random.seed(42)
    

    

    
    
    ####################################################################
    ###### read the discretized dataset
    ####################################################################
    
    file_name = '../dataset/artificial_matrix.csv'
    # file_name = '../dataset/artificial_matrix_worst.csv'
    
    # Load the discrete dataset
    df = pd.read_csv(file_name, index_col=0)

    
    print(df.head())
    
    
    # compute a UMAP from the data
    X = df.values.copy()
    # standard scaling before PCA
    X=StandardScaler().fit_transform(X)
    # PCA with ten principal components
    X_pca = PCA(n_components=10).fit_transform(X)

    reducer = umap.UMAP(min_dist=0.3,n_neighbors=50,spread=1.0)
    embedding = reducer.fit_transform(X_pca)
    
    fig, axs = plt.subplots(1,2)
    
    col = ['forestgreen' for _ in range(int(df.shape[0]/2))] + ['darkred' for _ in range(int(df.shape[0]/2))]
    
    axs[0].scatter(embedding[:, 0], embedding[:, 1], s=2, c=col)
    axs[0].set_aspect('equal', 'datalim')
    axs[0].set_title('projection UMAP')
    axs[0].set_xlabel('UMAP 1')
    axs[0].set_ylabel('UMAP 2')
    axs[0].legend(handles=[mpatches.Patch(color='forestgreen', label='exemples positifs'), mpatches.Patch(color='darkred', label='exemples négatifs')], loc='lower right')
    
    # PCA projection
    X = df.values.copy()
    X=StandardScaler().fit_transform(X)
    X_pca = PCA(n_components=2).fit_transform(X)
    axs[1].scatter(X_pca[:,0], X_pca[:,1], c=col, s=2)
    axs[1].set_title('projection PCA')
    axs[1].set_xlabel('PCA 1')
    axs[1].set_ylabel('PCA 2')
    axs[1].legend(handles=[mpatches.Patch(color='forestgreen', label='exemples positifs'), mpatches.Patch(color='darkred', label='exemples négatifs')])

    
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
    
    solver = Solver(inst)
    
    
    ############################################
    ### computation of supported front
    ############################################
    
    body_length = 20
    
    # create line equations from these points
    x1, y1 = 0.,0.
    x2, y2 = (inst.n_positives(), inst.n_negatives())
    a = 1.
    b = (x1-x2)/(y2-y1)
    c = -a*x1 - b*y1
    
    max_dist = a*0+b*y2+c
    
    # compute a line at distance 0.3
    p031 = [0,0]
    p031[1] = (0.3*max_dist-c)/b
    
    p032 = [0,y2]
    p032[0] = (0.3*max_dist-c-b*p032[1])/a
    
    # compute the list of distances to the diagonal
    atom_distance = [-1 for ind_atom in range(inst.n_atoms())]
    for ind_atom in range(inst.n_atoms()):
        score = inst.atom_score[ind_atom]
        atom_distance[ind_atom] = (a*score[0] + b*score[1] + c)/max_dist
    sorted_atoms = [ind for ind in range(inst.n_atoms())]
    sorted_atoms.sort(key=lambda ind:atom_distance[ind], reverse=True)
    
    
    stop_ind = 0
    for elt in sorted_atoms:
        if atom_distance[elt] < 0.3:
            break
        else:
            stop_ind += 1
    
    sorted_atoms = sorted_atoms[:stop_ind]
    
    # selection score = 0.75
    selection_score = 0.75
    
    selected_atoms, scores = solver.select_best_atoms_threshold(selection_score)

    
    print('n atoms selected: ', len(selected_atoms))

    
    body = selected_atoms
    
    # plot the body
    print([inst.get_atom(ind) for ind in body])
    
    nPos = inst.n_positives()
    nNeg = inst.n_negatives()
    
    ########################################################
    # compute pride rules
    pride_bodies = compute_pride_bodies(inst)
    pride_atoms = list(np.unique([ind_atom for pride_body in pride_bodies for ind_atom in pride_body]))
    print('n rules pride: ', len(pride_bodies))
    print('rule 0: ', pride_bodies[0])
    
    # inspection des indices des atomes des règles pride
    
    ###
    
    # number of atoms for pride rules
    # n_atoms = [len(elt) for elt in pride_bodies]
    # print('number of atoms pride: ', np.bincount(n_atoms))
    # print(sorted(n_atoms))
    
    # plot each atom score
    fig,ax = plt.subplots()
    
    ax.set_xlim((0,nPos))
    ax.set_ylim((0,nNeg))
    
    #ax.set_aspect('equal')
    
    ax.set_title('Atoms score')
    ax.set_xlabel('erreur positive')
    ax.set_ylabel('erreur négative')
    
    
    other_atoms = [ind_atom for ind_atom in range(inst.n_atoms()) if not ind_atom in body and not ind_atom in pride_atoms]
    
    # plot LOLH atoms
    ax.scatter([inst.atom_score[ind_atom][0] for ind_atom in body], [inst.atom_score[ind_atom][1] for ind_atom in body], marker='o', facecolor='darkorange', edgecolor='black', linewidth=0.5, zorder=2, label='atomes LOLH')
    
    # plot pride atoms
    ax.scatter([inst.atom_score[ind_atom][0] for ind_atom in pride_atoms], [inst.atom_score[ind_atom][1] for ind_atom in pride_atoms], marker='s', facecolor='dodgerblue', edgecolor='black', linewidth=0.5, zorder=2, label='atomes PRIDE')
    
    # other (not used) atoms
    ax.scatter([inst.atom_score[ind_atom][0] for ind_atom in range(inst.n_atoms())], [inst.atom_score[ind_atom][1] for ind_atom in range(inst.n_atoms())], alpha=0.5, marker='x', zorder=0, label='autres atomes')
    
    # score = 0.0
    ax.plot([0, nPos], [0, nNeg], '--', color='black', label='score = 0.0', lw=1.5)

    
    # compute the line corresponding to the selection score
    nPos = inst.n_positives()
    nNeg = inst.n_negatives()
    
    # line for score = 0.75
    A = (0,selection_score*nNeg)
    B = (nPos*(1.-selection_score), nNeg)
    ax.plot([A[0], B[0]], [A[1], B[1]], '--', color='forestgreen', label='score = 0.75', lw=1.5)
    
    ax.legend(loc='lower right')
    
    # compute a global score for each rule
    pride_rule_scores = []
    for pride_body in pride_bodies:
        score = (0, 0)
        for atom_index in pride_body:
            atom_score = inst.atom_score[atom_index]
            score = (score[0] + atom_score[0]/nPos, score[1] + atom_score[1]/nNeg)
        score = (score[0] / len(pride_body), score[1] / len(pride_body))
        pride_rule_scores.append(score)
    
    # compute the rule score for the LOLH rule
    lolh_rule_score = (0,0)
    for index in body:
        atom_score = inst.atom_score[index]
        lolh_rule_score = (lolh_rule_score[0] + atom_score[0]/nPos, lolh_rule_score[1] + atom_score[1]/nNeg)
    lolh_rule_score = (lolh_rule_score[0]/len(body), lolh_rule_score[1]/len(body))
    
    ############################################
    # plot the rule scores
    # display the normalized scores for all rules (between (0,0) and (n_positives, n_negatives))
    fig, ax = plt.subplots()
    
    ax.set_xlim((0, 1.))
    ax.set_ylim((0, 1.))
    ax.set_xlabel('score positif')
    ax.set_ylabel('score négatif')
    # ax.set_title('Comparaison des règles de PRIDE et LOLH')

    # diagonal indicating independant values
    ax.plot([0, 1], [0, 1], color='black', zorder=0) 
    
    # plot the pride rules
    pride_indexes = [ind for ind in range(len(pride_bodies))]
    pride_indexes.sort(key=lambda ind: len(pride_bodies[ind]))
    pride_rule_scores_sorted = [pride_rule_scores[ind] for ind in pride_indexes] # sort all the points
    longest_rule = np.max([len(body) for body in pride_bodies])
    colors = [len(pride_bodies[ind]) for ind in pride_indexes]
    cnorm = mcolors.Normalize(vmin=0, vmax=longest_rule)
    
    
    ax.scatter([score[0] for score in pride_rule_scores_sorted], [score[1] for score in pride_rule_scores_sorted], c=colors, norm=cnorm, cmap=plt.get_cmap('viridis'), marker='x',  zorder=1, label='règles PRIDE', facecolor='darkblue')
    cbar = fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap=plt.get_cmap('viridis')), ax=ax)
    cbar.set_label('longueur des règles PRIDE')
    ax.scatter(lolh_rule_score[0], lolh_rule_score[1], marker='*', s=70, color='chocolate', zorder=1, label='règle LOLH') # plot the optimized rule
    # plot the domination cone of the optimized rule
    ax.plot([lolh_rule_score[0], lolh_rule_score[0]], [lolh_rule_score[1], 0], linestyle='--', color='darkgrey', zorder=0) 
    ax.plot([lolh_rule_score[0], 1.], [lolh_rule_score[1], lolh_rule_score[1]], linestyle='--', color='darkgrey', zorder=0) 
    legend = ax.legend(loc='lower right')
    legend.legendHandles[0].set_color('darkcyan')
    
    ######################################################
    # create an histogram to visualize rule error over samples
    histo = histogram.Histogram(inst, body)
    # plot the histogram
    fig, ax = plt.subplots()
    #ax.set_title('Histograms for ' + instance_name + ' classification')
    visualizer.plot_histograms(ax, histo, True)
    
    ######################################################
    # histogramme de PRIDE à partir de l'apprentissage sur un sous jeu de données
    # histo = histogram.Histogram(inst, body)
    # # plot the histogram
    # fig, ax = plt.subplots()
    # #ax.set_title('Histograms for ' + instance_name + ' classification')
    # visualizer.plot_histograms(ax, histo, True)
    
    # histogramme de la meilleure règle PRIDE
    scores_pride = []
    ind_pride_bodies_sorted =  list(range(len(pride_bodies)))
    ind_pride_bodies_sorted.sort(key=lambda ind: pride_rule_scores[ind][1] - pride_rule_scores[ind][0], reverse=True)
    # print('pride scores sorted: ', [pride_rule_scores[ind][1] - pride_rule_scores[ind][0] for ind in ind_pride_bodies_sorted])

    ind_best_pride_body = ind_pride_bodies_sorted[0]
    print('best pride body: ', pride_bodies[ind_best_pride_body])
      
    # study atom indexes in pride bodies
    atom_indexes_pride = [inst.get_atom(ind) for body in pride_bodies for ind in body]
    print('pride body indexes (unique): ', np.unique(atom_indexes_pride, axis=0))

    
    # occ = np.bincount(atom_indexes_pride)
    # print('occurences: ')
    # for ind_atom in range(len(occ)):
    #     print(ind_atom, ', ', inst.atom_score[ind_atom][1]/nNeg - inst.atom_score[ind_atom][0]/nPos, ', ', occ[ind_atom])
    
    # fig, ax_ = plt.subplots()
    # ax_.scatter([inst.atom_score[ind_atom][0]/nPos for ind_atom in range(len(occ))], [inst.atom_score[ind_atom][1]/nPos for ind_atom in range(len(occ))], c=occ)
    # ax_.set_title('indexes des atomes de PRIDE')
    # ax_.set_xlabel('positive score')
    # ax_.set_ylabel('negative score')
    # cbar = fig.colorbar(cm.ScalarMappable(cmap=plt.get_cmap('viridis')), ax=ax_)
    # cbar.set_label('nb occurences')
    
    # histo = histogram.Histogram(inst, pride_bodies[ind_best_pride_body])
    # fig, ax = plt.subplots()
    # ax.set_title('Meilleure règle PRIDE')
    # visualizer.plot_histograms(ax, histo, True)
    
    
    
    ################################################################
    # histogramme de la plus longue règle PRIDE et de la meilleure règle PRIDE
    longest_pride = pride_bodies[0]
    for pride_body in pride_bodies:
        if len(pride_body) > len(longest_pride):
            longest_pride = pride_body
    histo1 = histogram.Histogram(inst, longest_pride)
    histo2 = histogram.Histogram(inst, pride_bodies[ind_best_pride_body])
    
    # plot the histogram
    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1,0.9]})
    
    fig.tight_layout(w_pad=2)
    
    # ax.set_title('Histogrammes pour PRIDE')
    visualizer.plot_histograms(axs[0], histo1, True)
    axs[0].set_title(r'Plus longue règle PRIDE')
    
    visualizer.plot_histograms(axs[1], histo2, True)
    axs[1].set_title(r'Meilleure règle PRIDE')
    
    
    
    ################################################################
    # visualisation des règles en multi-objectif
    body_length = 10
    biobj_score, biobj_bodies, biobj_weight = solver.compute_supported(body_length, 1, 1)
    
    biobj_score = [ (score[0]/float(inst.n_positives()), score[1]/float(inst.n_negatives())) for score in biobj_score]
    
    # sélection de deux règles "extrêmes"
    ind_left = 0
    for ind in range(len(biobj_score)):
        score = biobj_score[ind]
        if score[0] >= -0.2 and score[0] <= 0.2:
            if score[1] >= 6.4 and score[1] <= 6.6:
                print('left: ', ind)
                ind_left = ind
    print('left body: ', biobj_bodies[ind_left])
    left_score = biobj_score[ind_left]
    
    ind_right = 0
    for ind in range(len(biobj_score)):
        score = biobj_score[ind]
        if score[0] >= 2.9 and score[0] <= 3.0:
            if score[1] >= 9.8 and score[1] <= 10.0:
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
    
    ax.text(left_score[0] + 0.5, left_score[1]-0.2, r'$p_1$')
    ax.text(right_score[0] + 0.0, right_score[1]-0.6, r'$p_2$')
    
    delta = 0.5
    ax.set_xlim((0-delta, body_length+delta))
    ax.set_ylim((0-delta, body_length+delta))
    ax.set_aspect('equal')
    
    histo_left = histogram.Histogram(inst, biobj_bodies[ind_left])
    histo_right = histogram.Histogram(inst, biobj_bodies[ind_right])
    
    fig,axs = plt.subplots(1,2)
    fig.tight_layout(w_pad=3)
    axs[0].set_title(r'Histogramme de la règle $p_1$')
    visualizer.plot_histograms(axs[0], histo_left, True)
    axs[1].set_title(r'Histogramme de la règle $p_2$')
    visualizer.plot_histograms(axs[1], histo_right, True)
    
    
    ################################################################
    # display several rules with histograms
    
    fig, ax = plt.subplots()
    ax.set_xlim([0, nPos])
    ax.set_ylim([0, nNeg])
    ax.scatter([inst.atom_score[ind][0] for ind in range(inst.n_atoms())], [inst.atom_score[ind][1] for ind in range(inst.n_atoms())], marker='x', zorder=1, alpha=0.5)
    # ax.set_title('Visualisation de plusieurs règles logiques')
    # ax.plot([0, nPos], [0, nNeg], '--', color='grey', zorder=3, alpha=0.7)
    ax.set_xlabel('erreur positive')
    ax.set_ylabel('erreur négative')
    
    fig, axsh = plt.subplots(2,2)
    fig.tight_layout(h_pad=4)
    
    ##########################################
    def compute_A_B(score, nPos, nNeg):
        return (0, score*nNeg), ((1.-score)*nPos, nNeg)
    
    score1 = 0.3
    score2 = 0.4
    A1, B1 = compute_A_B(score1, nPos, nNeg)
    A2, B2 = compute_A_B(score2, nPos, nNeg)
    selected_atoms, scores = solver.select_best_atoms_threshold(0.0)
    atoms_sandwich = [selected_atoms[ind] for ind in range(len(selected_atoms)) if scores[ind] >= score1 and scores[ind] <= score2]
    np.random.shuffle(atoms_sandwich)
    nAtomsSel = 10
    body_worst = atoms_sandwich[:nAtomsSel]
    print('body: ', body_worst)
    ax.scatter([inst.atom_score[ind_atom][0] for ind_atom in body_worst], [inst.atom_score[ind_atom][1] for ind_atom in body_worst], marker='x', color='orangered', zorder=4)
    # histogramme
    histo = histogram.Histogram(inst, body_worst)
    axsh[0,0].set_title('Règle sous-optimale')
    visualizer.plot_histograms(axsh[0,0], histo, True)
    
    
    ax.fill([A1[0], B1[0], B2[0], A2[0]], [A1[1], B1[1], B2[1], A2[1]], zorder=3, color='green', alpha=0.3, label='règle sous-optimale')
    ax.plot([A1[0], B1[0]], [A1[1], B1[1]], '--', color='k', zorder=3)
    ax.plot([A2[0], B2[0]], [A2[1], B2[1]], '--', color='k', zorder=3)
    
    score1 = -0.8
    score2 = -0.7
    A1, B1 = compute_A_B(score1, nPos, nNeg)
    
    A2, B2 = compute_A_B(score2, nPos, nNeg)
    
    selected_atoms, scores = solver.select_best_atoms_threshold(-1.)
    atoms_sandwich = [selected_atoms[ind] for ind in range(len(selected_atoms)) if scores[ind] >= score1 and scores[ind] <= score2]
    np.random.shuffle(atoms_sandwich)
    nAtomsSel = 10
    body_worst = atoms_sandwich[:nAtomsSel]
    print('body: ', body_worst)
    ax.scatter([inst.atom_score[ind_atom][0] for ind_atom in body_worst], [inst.atom_score[ind_atom][1] for ind_atom in body_worst], marker='x', color='orangered', zorder=4)
    # histogramme
    histo = histogram.Histogram(inst, body_worst)
    axsh[0,1].set_title('Règle inverse')
    visualizer.plot_histograms(axsh[0,1], histo, True)
    
    ax.fill([A1[0], B1[0], B2[0], A2[0]], [A1[1], B1[1], B2[1], A2[1]], zorder=3, color='orange', alpha=0.3, label='règle inverse')
    ax.plot([A1[0], B1[0]], [A1[1], B1[1]], '--', color='k', zorder=3)
    ax.plot([A2[0], B2[0]], [A2[1], B2[1]], '--', color='k', zorder=3)
    
    score1 = -0.05
    score2 = 0.05
    A1, B1 = compute_A_B(score1, nPos, nNeg)
    A2, B2 = compute_A_B(score2, nPos, nNeg)
    selected_atoms, scores = solver.select_best_atoms_threshold(-1.)
    atoms_sandwich = [selected_atoms[ind] for ind in range(len(selected_atoms)) if scores[ind] >= score1 and scores[ind] <= score2]
    np.random.shuffle(atoms_sandwich)
    nAtomsSel = 10
    body_worst = atoms_sandwich[:nAtomsSel]
    print('body: ', body_worst)
    ax.scatter([inst.atom_score[ind_atom][0] for ind_atom in body_worst], [inst.atom_score[ind_atom][1] for ind_atom in body_worst], marker='x', color='orangered', zorder=4, label='atomes sélectionnés')
    # histogramme
    histo = histogram.Histogram(inst, body_worst)
    axsh[1,0].set_title('Règle inappropriée')
    visualizer.plot_histograms(axsh[1,0], histo, True)
    
    
    ax.fill([A1[0], B1[0], B2[0], A2[0]], [A1[1], B1[1], B2[1], A2[1]], zorder=3, color='red', alpha=0.3, label='règle inappropriée')
    ax.plot([A1[0], B1[0]], [A1[1], B1[1]], '--', color='k', zorder=3)
    ax.plot([A2[0], B2[0]], [A2[1], B2[1]], '--', color='k', zorder=3)
    
    atomes = [ind for ind in range(inst.n_atoms())]
    np.random.shuffle(atomes)
    body_random = atomes[:10]
    histo = histogram.Histogram(inst, body_random)
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
    
    
    ################################################################
    # K-fold validation on PRIDE and LOLH
    k = 10
    
    # create nsamples/k subset of the data
    pos_indexes = inst._pos_samples
    pos_sub_indexes = [[] for _ in range(k)]
    for ind in range(len(pos_indexes)):
        k_prime = int((ind/nPos)*k)
        pos_sub_indexes[k_prime].append(pos_indexes[ind])
        
    neg_indexes = inst._neg_samples
    neg_sub_indexes = [[] for _ in range(k)]
    for ind in range(len(neg_indexes)):
        k_prime = int((ind/nNeg)*k)
        neg_sub_indexes[k_prime].append(neg_indexes[ind])
    
    LOLH_perf = [[],[]]
    PRIDE_perf = [[],[]]
    
    for ind_k in range(1):
    
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
        LOLH_body, scores = solver.select_best_atoms_threshold(0.75)
        LOLH_body = [inst_train.get_atom(atom_index) for atom_index in LOLH_body]
        LOLH_body = [inst_test.get_atom_index(atom) for atom in LOLH_body]
        
        # compute PRIDE bodies
        PRIDE_bodies = compute_pride_bodies(inst_train)
        PRIDE_bodies = [[inst_train.get_atom(atom_index) for atom_index in pride_body] for pride_body in PRIDE_bodies]
        PRIDE_bodies = [[inst_test.get_atom_index(atom) for atom in PRIDE_body] for PRIDE_body in PRIDE_bodies]
        
        
        # atoms in the train and test instance
        atoms_train_test = [inst_train.get_atom(atom_index) for atom_index in range(inst_train.n_atoms()) if inst_test.has_atom(inst_train.get_atom(atom_index))]
        col = []
        for atom in atoms_train_test:
            ind_atom_train = inst_train.get_atom_index(atom)
            score = -inst_train.atom_score[ind_atom_train][0]*len(train_neg)+inst_train.atom_score[ind_atom_train][1]*len(train_pos)
            col.append(score)
                
        # fig, ax = plt.subplots()
        # ax.scatter([inst_test.atom_score[inst_test.get_atom_index(atom)][0] for atom in atoms_train_test], [inst_test.atom_score[inst_test.get_atom_index(atom)][1] for atom in atoms_train_test], marker='x', c=col)
        # ax.set_xlim([0, len(test_pos)])
        # ax.set_ylim([0, len(test_neg)])
        # ax.plot([0, len(test_pos)], [0, len(test_neg)], '--', color='k')
    
        # histogramme de LOLH
        # histo = histogram.Histogram(inst_test, LOLH_body)
        # fig, ax = plt.subplots()
        # ax.set_title('histo LOLH données test')
        # visualizer.plot_histograms(ax, histo, True)
        
        # LOLH performance
        pos_match = 0
        neg_match = 0
        for pos in test_pos:
            matching_error = 0
            for ind_atom in LOLH_body:
                atom = inst_test.get_atom(ind_atom)
                if df_test[atom[0]][pos] != atom[1]:
                    matching_error += 1
            if matching_error < len(LOLH_body)/2.:
                pos_match += 1
        for neg in test_neg:
            matching_error = 0
            for ind_atom in LOLH_body:
                atom = inst_test.get_atom(ind_atom)
                if df_test[atom[0]][neg] != atom[1]:
                    matching_error += 1
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
            match = False
            for pride_body in PRIDE_bodies:
                matching_error = 0
                for ind_atom in pride_body:
                    atom = inst_test.get_atom(ind_atom)
                    if df_test[atom[0]][pos] != atom[1]:
                        matching_error += 1
                if matching_error == 0:
                    match = True
                    pos_match += 1
                    break
        for neg in test_neg:
            match = False
            for pride_body in PRIDE_bodies:
                matching_error = 0
                for ind_atom in pride_body:
                    atom = inst_test.get_atom(ind_atom)
                    if df_test[atom[0]][neg] != atom[1]:
                        matching_error += 1
                if matching_error == 0:
                    match = True
                    neg_match += 1
                    break
        PRIDE_perf[0].append(pos_match/len(test_pos))
        PRIDE_perf[1].append(neg_match/len(test_neg))
        
    print(LOLH_perf)
    print(PRIDE_perf)
    
    return

# singlecell_matrix_experiment()

# artificial_data_experiment()

LOLH_artificial_data()
