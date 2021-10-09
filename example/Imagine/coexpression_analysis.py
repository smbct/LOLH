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


#------------------------------------------------------------------------------
#                Load the gene clusters from text files
#------------------------------------------------------------------------------
def load_clusters(filename):
    file = open(filename, 'r')
    content = file.read().splitlines()
    file.close()
    gene_clusters = []
    for line in content:
        gene_cluster = []
        atoms_str = line.split(' ')
        for elt in atoms_str:
            tokens = elt.split('_')
            if len(tokens) == 2:
                atom = (tokens[0], int(tokens[1]))
                gene_cluster.append(atom)
        gene_clusters.append(gene_cluster)
    return gene_clusters




#------------------------------------------------------------------------------
#                Print the gene clusters formatted in LaTex
#------------------------------------------------------------------------------
def print_cluster_latex(instance, cluster):

    # output the clusters formatted in LaTex
    res = ''
    for elt in cluster:
        atom_index = instance.get_atom_index(elt)
        val = elt[1]
        nval = instance.n_values[elt[0]]
        color_str = None
        if nval == 2:
            if val == 0:
                color_str = 'BrickRed'
            elif val == 1:
                color_str = 'OliveGreen'
        elif nval == 3:
            if val == 0:
                color_str = 'BrickRed'
            elif val == 1:
                color_str = 'Goldenrod'
            elif val == 2:
                color_str = 'OliveGreen'
        elif nval == 4:
            if val == 0:
                color_str = 'BrickRed'
            elif val == 1:
                color_str = 'Orange'
            elif val == 2:
                color_str = 'SpringGreen'
            elif val == 3:
                color_str = 'OliveGreen'
        res += '\\textcolor{' + color_str + '}{$\\boldsymbol{\\textit{' + elt[0] + '}' + '_{' + str(elt[1]) + '/' + str(nval-1) + '}}$}, '
    print(res)

    return



#------------------------------------------------------------------------------
#                Print the gene clusters
#------------------------------------------------------------------------------
def print_gene_clusters(gene_clusters):
    for ind_cluster in range(len(gene_clusters)):
        print('gene cluster ', ind_cluster, ': ')
        str_clusters = ''
        for atom in gene_clusters[ind_cluster]:
            str_clusters += atom[0] + '_' + str(atom[1]) + ', '
        print(str_clusters)
        print('\n\n')
    return

#------------------------------------------------------------------------------
#                Plot the rule matching error for a body
#------------------------------------------------------------------------------
def plot_rule_error(df_discrete, df_coordinates, instance, body, ax):

    histo = histogram.Histogram(instance, body, histogram.Histogram_type.GLOBAL)

    plasma = cm.get_cmap('plasma_r', 256)

    cnorm = mcolors.Normalize(vmin=0, vmax=len(histo.positive_histogram))

    cell_score = {barcode : error for error in range(len(histo.positive_histogram)) for barcode in histo.positive_histogram[error]}

    # sort all the cells according to the score
    barcodes = [index for index in df_coordinates.index]
    barcodes.sort(key = lambda elt: cell_score[elt], reverse=True)
    df_coordinates_sorted = df_coordinates.loc[barcodes]

    col = [cell_score[barcode] for barcode in df_coordinates_sorted.index]
    ax.scatter(df_coordinates_sorted['UMAP_1'].values, df_coordinates_sorted['UMAP_2'].values, c=col, cmap=plasma, norm=cnorm, s=3)

    # remove UMAP coordinates
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
    cbar.set_label('Matching error')

    return

#------------------------------------------------------------------------------
#                Plot the cell scores for each cluster
#------------------------------------------------------------------------------
def plot_cell_scores(df_discrete, df_coordinates, gene_clusters, cluster_selection = None, n_col = 3, title=None):

    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)

    if cluster_selection == None:
        cluster_indexes = [ind for ind in range(len(gene_clusters))]
    else:
        cluster_indexes =  cluster_selection

    n_line = math.ceil(len(cluster_indexes)/n_col)
    fig, axs = plt.subplots(n_line, n_col)

    ind_plot = 0
    for ind_line in range(n_line):
        for ind_col in range(n_col):
            if ind_plot < len(cluster_indexes):
                ax = axs[ind_line, ind_col]
                ind_cluster = cluster_indexes[ind_plot]
                body = [instance.get_atom_index(atom) for atom in gene_clusters[ind_cluster] ]
                plot_rule_error(df_discrete, df_coordinates, instance, body, ax)
                ax.set_title('cluster ' + str(ind_cluster))
                ind_plot += 1

    if title != None:
        fig.suptitle(title)

    return

#------------------------------------------------------------------------------
#                create a dataframe with the gene clusters identity
#------------------------------------------------------------------------------
def create_dataframe_clusters(filename, gene_clusters):

    genes = np.concatenate(gene_clusters)
    data = [ [int(elt[1]) ] for elt in genes]
    for ind_gene in range(len(genes)):
        for ind_cluster in range(len(gene_clusters)):
            gene_cluster = gene_clusters[ind_cluster]
            if (genes[ind_gene][0], int(genes[ind_gene][1])) in gene_cluster:
                data[ind_gene].append(1)
            else:
                data[ind_gene].append(0)
    df_gene_clusters = pd.DataFrame(data, index=[elt[0] for elt in genes], columns=['gene_value'] + ['cluster_'+str(ind) for ind in range(len(gene_clusters))])
    df_gene_clusters.to_csv(filename)

    return




#------------------------------------------------------------------------------
#                create a dataframe containing the matching error for all cells on the clusters
#------------------------------------------------------------------------------
def create_cell_score_dataframe(filename, df_discrete, gene_clusters):

    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)

    data = [[] for _ in df_discrete.index]
    barcode_to_ind = { df_discrete.index[ind]:ind for ind in range(len(df_discrete.index))}
    for ind_cluster in range(len(gene_clusters)):
        print('gene cluster index: ', ind_cluster)
        body = [instance.get_atom_index(atom) for atom in gene_clusters[ind_cluster] ]
        histo = histogram.Histogram(instance, body, histogram.Histogram_type.GLOBAL)
        for error in range(len(body)+1):
            cell_score = 1.-error/len(body)
            for barcode in histo.positive_histogram[error]:
                data[barcode_to_ind[barcode]].append(cell_score)
    df_cell_score = pd.DataFrame(data, index=df_discrete.index, columns=['gene_cluster_' + str(ind) for ind in range(len(gene_clusters))])

    df_cell_score.to_csv(filename)

    return



#------------------------------------------------------------------------------
#                create a discrete dataframe with a subset of the cells
#------------------------------------------------------------------------------
def create_sub_dataframe(filename, df_discrete, selected_cells):

    df_sub_discrete = df_discrete.copy()
    df_sub_discrete = df_sub_discrete.loc[selected_cells]
    df_sub_discrete = df_sub_discrete.loc[:, (df_sub_discrete != df_sub_discrete.iloc[0]).any()] # remove potential constant columns

    print(df_sub_discrete)

    # print('filename: ', filename)

    df_sub_discrete.to_csv(filename)

    return

#------------------------------------------------------------------------------
#                compose the cells into two subset based on the cluster matching error
#------------------------------------------------------------------------------
def decompose_cluster(df_coordinates, instance, gene_clusters, cluster_index, error_threshold, all_cells = False):

    # plot global histogram from cluster 6 and select cells with a small matching error
    body = [instance.get_atom_index(atom) for atom in gene_clusters[cluster_index] ]

    histo = histogram.Histogram(instance, body, histogram.Histogram_type.GLOBAL)
    fig,ax = plt.subplots()
    visualizer.plot_global_histogram(ax, histo)
    ax.set_ylabel('cell proportion')
    title = 'Cell matching error for gene cluster ' + str(cluster_index)
    if all_cells:
        title += ' (on all cells)'
    ax.set_title(title)
    selected_cells = [cell for error in range(error_threshold+1) for cell in histo.positive_histogram[error]]
    selected_cells
    remaining_cells = [barcode for barcode in df_coordinates.index if not barcode in selected_cells]

    # plot the selected cells from cluster 6 in the UMAP
    fig, ax = plt.subplots()
    ax.scatter(df_coordinates.loc[remaining_cells]['UMAP_1'].values, df_coordinates.loc[remaining_cells]['UMAP_2'].values, s=3, c='black')
    ax.scatter(df_coordinates.loc[selected_cells]['UMAP_1'].values, df_coordinates.loc[selected_cells]['UMAP_2'].values, s=3, c='red')

    title = 'Cells selected from matching error on cluster ' + str(cluster_index)
    if all_cells:
        title += ' (on all cells)'
    ax.set_title(title)

    # export the cells corresponding to cluster 6 in a csv file
    df_selection = pd.DataFrame([1 for _ in selected_cells] + [0 for _ in remaining_cells], index=selected_cells+remaining_cells, columns=['selection'])
    filename = '../../coexpression/myeloid_c' + str(cluster_index) + '_selection'
    if all_cells:
        filename += '_all_cells'
    filename += '.csv'
    df_selection.to_csv(filename)

    return

#------------------------------------------------------------------------------
#                process the global clusters
#------------------------------------------------------------------------------
def process_global_clusters():

    # read the discrete matrix
    filename = '../../dataset/discrete_dataset.csv'
    df_discrete = pd.read_csv(filename, index_col=0)
    print(df_discrete)

    # read the embedding coordinates
    df_coordinates = pd.read_csv('../../dataset/UMAP_coordinates.csv', index_col=0)


    # read the gene clusters
    gene_clusters = load_clusters('../../coexpression/gene_clusters.txt')

    print('n clusters: ', len(gene_clusters))

    # create an artifial instance
    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)

    # output the clusters formatted in LaTex
    # for cluster_index in range(len(gene_clusters)):
    #     cluster = gene_clusters[cluster_index]
    #     print('gene cluster ', cluster_index, '(', len(cluster), ' atoms)')
    #     print_cluster_latex(instance, gene_clusters[cluster_index])
    #     print('\n\n')

    # Display the cell score for each gene cluster on the UMAP
    plot_cell_scores(df_discrete, df_coordinates, gene_clusters, None, 3, None)


    # display gene cluster 4 cell score in a histogram
    ind_cluster = 4
    body = [instance.get_atom_index(atom) for atom in gene_clusters[ind_cluster] ]
    histo = histogram.Histogram(instance, body, histogram.Histogram_type.GLOBAL)


    # visualization of the matching error histogram of cluster #4 on all the cells
    fig,ax = plt.subplots()
    values = [ error for error in range(len(histo.positive_histogram)) for _ in histo.positive_histogram[error] ]
    ax.hist(values, 50, density=True, edgecolor='black')
    ax.set_ylabel('cell proportion')
    ax.set_xlabel('matching error')
    ax.set_title('Cell matching error for gene cluster ' + str(ind_cluster))


    # select the cells that have the lower score
    threshold = 193
    selected_cells = []
    for error in range(threshold+1):
        selected_cells += histo.positive_histogram[error]

    # plot a UMAP with the selected cells
    fig, ax = plt.subplots()
    ax.set_xlabel('UMAP_1')
    ax.set_ylabel('UMAP_2')
    ax.set_title('Selected cells from gene cluster ' + str(ind_cluster))
    colors = ['red' if barcode in selected_cells else 'black' for barcode in df_coordinates.index]
    ax.scatter(df_coordinates['UMAP_1'].values, df_coordinates['UMAP_2'].values, c=colors, s=3)
    ax.set_title('Cells selected through the matching error from cluster ' + str(ind_cluster))

    # print the gene clusters
    print_gene_clusters(gene_clusters)

    # creation of a DataFrame with the gene clusters
    #create_dataframe_clusters('../../coexpression/gene_clusters.csv', gene_clusters)

    # creation of a dataframe with the cell scores (between 0 and 1)
    #create_cell_score_dataframe('../../coexpression/gene_cluster_cell_scores.csv', df_discrete, gene_clusters)

    # create a sub dataset with the selected cells (presumably Myeloids)
    create_sub_dataframe('../../dataset/sub_dataset_discrete.csv', df_discrete, selected_cells)

    plt.show()

    return


#------------------------------------------------------------------------------
#                process the sub network clusters
#------------------------------------------------------------------------------
def process_sub_network_clusters():

    # read the sub discrete matrix (with only cells selected from the first network)
    filename = '../../dataset/sub_dataset_discrete.csv'
    df_discrete_sub = pd.read_csv(filename, index_col=0)
    print(df_discrete_sub.head())

    # read the complete discrete matrix (to visualize the cluster matching error on all the cells)
    filename = '../../dataset/discrete_dataset.csv'
    df_discrete = pd.read_csv(filename, index_col=0)

    # read the embedding coordinates
    df_coordinates = pd.read_csv('../../dataset/UMAP_coordinates.csv', index_col=0)
    df_coordinates_sub = df_coordinates.loc[df_discrete_sub.index] # select only cells in the sub dataset

    # read the gene clusters
    gene_clusters = load_clusters('../../coexpression/sub_network_gene_clusters.txt')


    ###########################################################################
    # Manual sub selection of the Myeloids
    selected_indexes = []
    for index in df_coordinates_sub.index:
        pos = (df_coordinates_sub['UMAP_1'][index], df_coordinates_sub['UMAP_2'][index])
        if pos[0] >= -8 and pos[0] <= 8.0:
            if pos[1] >= -16 and pos[1] <= -11:
                selected_indexes.append(index)
    df_coordinates_sub_manual = df_coordinates_sub.loc[selected_indexes]
    ###########################################################################

    # create a artificial instance to compute the matching error
    instance = Instance.create_random_instance(df_discrete_sub.copy(deep=False), 0.5)

    selected_clusters = [1,4,5,8,9,10]

    # output the clusters formatted in LaTex
    # for cluster_index in selected_clusters:
    #     cluster = gene_clusters[cluster_index]
    #     print('gene cluster ', cluster_index, '(', len(cluster), ' atoms)')
    #     print_cluster_latex(instance, gene_clusters[cluster_index])
    #     print('\n\n')

    # Display the cell score for each gene cluster on the UMAP
    plot_cell_scores(df_discrete_sub, df_coordinates_sub_manual, gene_clusters, selected_clusters, 3, 'Clusters matching error on the myeloids')

    # Same visualization on all the cells this time
    plot_cell_scores(df_discrete, df_coordinates, gene_clusters, selected_clusters, 3, 'Clusters matching error on all the cells')

    # analysis of cluster 8 matching error
    # cluster_index = 8
    # error_threshold = 6
    # decompose_cluster(df_coordinates_sub, instance, gene_clusters, cluster_index, error_threshold, False)

    instance_global = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5) # create a artificial instance on al cells

    # analysis of cluster 8 matching error on all cells
    # cluster_index = 8
    # error_threshold = 8
    # decompose_cluster(df_coordinates, instance_global, gene_clusters, cluster_index, error_threshold, True)


    # analyss of cluster 8: selection of cells and decomposition into two groups: cells in the sub dataset and other cells
    cluster_index = 8
    error_threshold = 8
    body = [instance_global.get_atom_index(atom) for atom in gene_clusters[cluster_index] ]
    histo = histogram.Histogram(instance_global, body, histogram.Histogram_type.GLOBAL)
    selected_cells = [cell for error in range(error_threshold+1) for cell in histo.positive_histogram[error]]
    other_cells = [barcode for barcode in df_discrete.index if not barcode in selected_cells]
    positive_cells = [barcode for barcode in selected_cells if barcode in df_discrete_sub.index]
    negative_cells = [barcode for barcode in selected_cells if not barcode in positive_cells]
    fig,ax = plt.subplots()
    ax.scatter(df_coordinates.loc[other_cells]['UMAP_1'].values, df_coordinates.loc[other_cells]['UMAP_2'].values, c='k', s=3)
    ax.scatter(df_coordinates.loc[positive_cells]['UMAP_1'].values, df_coordinates.loc[positive_cells]['UMAP_2'].values, c='r', s=3)
    ax.scatter(df_coordinates.loc[negative_cells]['UMAP_1'].values, df_coordinates.loc[negative_cells]['UMAP_2'].values, c='b', s=3)
    ax.set_title('Selected cells from cluster 8: myeloids (red) vs other cells (blue)')
    ax.set_xlabel('UMAP_1')
    ax.set_ylabel('UMAP_2')
    df_selection = pd.DataFrame([1 for _ in positive_cells] + [0 for _ in negative_cells], index=positive_cells+negative_cells, columns=['selection'])
    filename = '../../coexpression/myeloid_c8_selection_myeloid_vs_others.csv'
    df_selection.to_csv(filename)

    # analysis of cluster 10 matching error
    cluster_index = 10
    error_threshold = 38
    decompose_cluster(df_coordinates_sub, instance, gene_clusters, cluster_index, error_threshold, False)

    # print the gene clusters
    print_gene_clusters(gene_clusters)

    # creation of a DataFrame with the gene clusters
    create_dataframe_clusters('../../coexpression/sub_network_gene_clusters.csv', gene_clusters)

    # creation of a dataframe with the cell scores (between 0 and 1)
    create_cell_score_dataframe('../../coexpression/sub_network_gene_cluster_cell_scores.csv', df_discrete, gene_clusters)

    plt.show()

    return

# process_global_clusters()

# process_sub_network_clusters()
