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

from network import GData
from network import Graph

#-------------------------------------------------------------------------------
def load_graph_file(filename):
    file = open(filename, 'r')
    graph = file.read().splitlines()
    file.close()
    return graph

#-------------------------------------------------------------------------------
def save_graph_file(filename, graph_reduced):
    file = open(filename, 'w')
    for ind_line in range(len(graph_reduced)):
        line = graph_reduced[ind_line]
        for ind in range(len(line)):
            file.write(line[ind])
            if ind < len(line)-1:
                file.write(' ')
        file.write('\n')
    file.close()
    return

#-------------------------------------------------------------------------------
def reduce_graphe(graph, threshold):

    graph_bis = []
    for line in graph:
        tokens = line.split(' ')
        graph_bis.append([tokens[i] for i in range(4)])
        ind = 4
        while ind < len(tokens):
            val = float(tokens[ind+2])
            if val >= threshold:
                graph_bis[-1] += [tokens[ind+i] for i in range(5)]
            ind += 5
    return graph_bis


#-------------------------------------------------------------------------------
def create_regulatory_graph():

    filename = '../../dataset/Imagine/regulatory_network03.txt'
    graph = load_graph_file(filename)

    # create a new file and truncate some edges
    threshold = 0.5
    graph_reduced = reduce_graphe(graph, threshold)

    filename = '../../dataset/Imagine/regulatory_network05.txt'
    save_graph_file(filename, graph_reduced)

    filename = '../../dataset/Imagine/regulatory_network05.txt'
    ncell_min = 50
    score_min = 0.5
    data = GData()
    data.load_from_file(filename, ncell_min, score_min)

    # construction
    print('Construction of the regulatory graph')
    exclude_mt_rp = True
    filter_edges = False
    graph = Graph('regulatory network')
    graph.create_from_gdata_regul(data, exclude_mt_rp, filter_edges = False)

    # compute vertex 2d position
    print('compute 2d positions for the vertices')
    graph.compute_positions()

    # clustering
    print('clustering of the graph')
    res = 1
    clusters = graph.compute_clusters(res)

    return graph


#-------------------------------------------------------------------------------
def plot_cluster_umap(df_discrete, df_coordinates, instance, body, ax):

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

    # squared plots
    # ax.set_aspect((ax.get_ylim()[1]-ax.get_ylim()[0])/(ax.get_xlim()[1]-ax.get_xlim()[0]))
    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))

    cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
    # cbar.set_label('Matching error')
    cbar.set_label('Erreur de couv.')

    return


#-------------------------------------------------------------------------------
def plot_gene_umap(df_normalized, df_coordinates, gene, ax):

    plasma = cm.get_cmap('plasma', 256)

    cell_values_dic = {barcode : df_normalized[gene][barcode] for barcode in df_normalized.index}
    cell_values = list(cell_values_dic.values())

    barcodes = [index for index in df_coordinates.index]
    barcodes.sort(key = lambda elt: cell_values_dic[elt], reverse=False)
    df_coordinates_sorted = df_coordinates.loc[barcodes]

    cnorm = mcolors.Normalize(vmin=min(cell_values), vmax=max(cell_values))

    col = [cell_values_dic[barcode] for barcode in df_coordinates_sorted.index]

    ax.scatter(df_coordinates_sorted['UMAP_1'].values, df_coordinates_sorted['UMAP_2'].values, c=col, cmap=plasma, norm=cnorm, s=3)

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
    cbar.set_label('gene value')

    return



#-------------------------------------------------------------------------------
def global_analysis():


    # create_regulatory_graph()

    # load the graph from an existing file
    filename = '../../dataset/Imagine/regulatory_network_processed.txt'
    graph = Graph('regulatory network')
    graph.load_from_file(filename)


    # load the discrete dataset
    filename = '../../dataset/Imagine/discrete_matrix.csv'
    df_discrete = pd.read_csv(filename, index_col = 0)
    print(df_discrete.head())

    # load the cell types
    filename = '../../dataset/Imagine/cell_types.csv'
    df_celltypes = pd.read_csv(filename, index_col = 0)
    print(df_celltypes.head())

    # load the UMAP representation
    filename = '../../dataset/Imagine/umap_coordinates.csv'
    df_umap = pd.read_csv(filename, index_col = 0)
    print(df_umap.head())

    # load the normalized dataset
    filename = '../../dataset/Imagine/normalized_matrix.csv'
    df_normalized = pd.read_csv(filename, index_col = 0)
    df_normalized = df_normalized.T
    print(df_normalized.head())

    # display
    print('display the (directed) graph')
    col_option = '01_colors'
    arrows = True
    cluster_size_limit = 20
    fig, ax = plt.subplots()
    graph.plot(ax, col_option, arrows, cluster_size_limit)

    col_option = 'clustering_colors'
    fig, ax = plt.subplots()
    graph.plot(ax, col_option, arrows, cluster_size_limit)


    ################################################################################
    #### load and process the coexpression network
    # input_file = '../../dataset/Imagine/coexpression/coexpression_network.txt'
    # n_cell_min = 200
    # score_min = 0.35
    # louvain_param = 0.7
    # output_file = '../../dataset/Imagine/coexpression/coexpression_graph.txt'
    # codata = GData()
    # codata.load_from_file(input_file, n_cell_min, score_min)
    # cograph = Graph('coexpression network')
    # # build graph from raw data, exclude mitochondrial and ribosomal genes
    # exclude_mt_rp = True
    # filter_edges = True
    # cograph.create_from_gdata(codata, exclude_mt_rp, filter_edges)
    # cograph.compute_clusters(louvain_param)
    # cograph.compute_positions()
    # cograph.save(output_file)

    cograph = Graph('coexpression network')
    cograph.load_from_file('../../dataset/Imagine/coexpression/coexpression_graph.txt')
    fig, ax = plt.subplots()
    ax.set_title('coexpression graph')
    cograph.plot(ax, 'clustering_colors', False, 10)

    # save the dynamical graph
    # filename = '../../dataset/Imagine/regulatory_network_processed.txt'
    # graph.save(filename)

    clusters = {}
    for cluster in graph.clusters:
        if len(graph.clusters[cluster]) >= 20:
            clusters[cluster] = graph.clusters[cluster]
    # print(clusters)

    print('plot the clusters')
    # create a fake instance
    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    # cluster_index = list(clusters.keys())[0]

    selected_clusters = [6,7,8,9,11,20]

    ncol = 2
    nrows = int(len(selected_clusters)/float(ncol))
    fix, axs = plt.subplots(nrows, ncol)
    axs = axs.flat
    ind_plot = 0
    for cluster_index in selected_clusters:
        ax = axs[ind_plot]
        body =  [ instance.get_atom_index(graph.atoms[index]) for index in clusters[ cluster_index ] ]
        # fig, ax = plt.subplots()
        ax.set_title('cluster ' + str(cluster_index))
        plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
        ind_plot += 1


    # plot a single gene
    gene = 'PPBP'
    gene = 'S100A4'
    gene = 'ANXA1'
    gene = 'DERL3'

    fig, ax = plt.subplots()
    ax.set_title('gene ' + gene)
    plot_gene_umap(df_normalized, df_umap, gene, ax)

    plt.show()

    return

#-------------------------------------------------------------------------------
def sub_analysis():

    # create the sub dataset limited to the t cell
    # load the discrete dataset
    # filename = '../../dataset/Imagine/discrete_matrix.csv'
    # df_discrete = pd.read_csv(filename, index_col = 0)
    # # load the cell types
    # filename = '../../dataset/Imagine/cell_types_macro.csv'
    # df_celltypes = pd.read_csv(filename, index_col = 0)
    # # create a sub dataset
    # T_cells = df_celltypes.loc[df_celltypes['cellType_macro']=='T'].index
    # sub_dataset = df_discrete.loc[T_cells]
    # sub_dataset.to_csv('../../dataset/Imagine/discrete_matrix_T_cells.csv')

    # create the subset of transitions limited to the t cells
    # transitions_df = pd.read_csv('../../dataset/Imagine/transitions.csv', index_col=0)
    # print(transitions_df.head())
    # filename = '../../dataset/Imagine/discrete_matrix_T_cells.csv'
    # df_T = pd.read_csv(filename, index_col = 0)
    # T_cells = list(df_T.index)
    # selected_transitions = []
    # for ind in transitions_df.index:
    #     if transitions_df['T-1'][ind] in T_cells and transitions_df['T'][ind] in T_cells:
    #         selected_transitions.append(ind)
    # transitions_df_sub = transitions_df.loc[selected_transitions]
    # transitions_df_sub.to_csv('../../dataset/Imagine/transitions_T.csv')

    filename = '../../dataset/Imagine/T_regulatory_network.txt'

    n_cell_min = 20
    score_min = 0.35
    louvain_param = 0.7
    data = GData()
    data.load_from_file(filename, n_cell_min, score_min)
    graph = Graph('regulatory network (T cells)')
    # build graph from raw data, exclude mitochondrial and ribosomal genes
    exclude_mt_rp = True
    filter_edges = True
    graph.create_from_gdata(data, exclude_mt_rp, filter_edges)
    graph.compute_clusters(louvain_param)
    graph.compute_positions()

    fig, ax = plt.subplots()
    ax.set_title('T cels regulatory graph')
    graph.plot(ax, 'clustering_colors', True, 10)

    fig, ax = plt.subplots()
    ax.set_title('T cels regulatory graph')
    graph.plot(ax, '01_colors', True, 10)

    # graph.save(output_file)

    plt.show()

    return


# global_analysis()

sub_analysis()
