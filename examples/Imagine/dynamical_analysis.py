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

    filename = '../../dataset/Imagine/dynamics/regulatory_network03.txt'
    graph = load_graph_file(filename)

    # create a new file and truncate some edges
    threshold = 0.5
    graph_reduced = reduce_graphe(graph, threshold)

    filename = '../../dataset/Imagine/dynamics/regulatory_network05.txt'
    save_graph_file(filename, graph_reduced)

    filename = '../../dataset/Imagine/dynamics/regulatory_network05.txt'
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
    ax.scatter(df_coordinates_sorted['UMAP_1'].values, df_coordinates_sorted['UMAP_2'].values, c=col, cmap=plasma, norm=cnorm, s=1)

    # remove UMAP coordinates
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    # squared plots
    # ax.set_aspect((ax.get_ylim()[1]-ax.get_ylim()[0])/(ax.get_xlim()[1]-ax.get_xlim()[0]))
    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))

    cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
    cbar.set_label('Matching error')
    # cbar.set_label('Erreur de couv.')

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

#------------------------------------------------------------------------------
#                Print the gene clusters formatted in LaTex
#------------------------------------------------------------------------------
def print_cluster_latex(instance, cluster):

    print(cluster)

    cluster_values = { (0,3):0, (0,2):1, (0,1):2, (1,3):3, (1,2):4, (2,3):5, (1,1):6, (2,2):7, (3,3):9}

    ordered_cluster = cluster.copy()
    ordered_cluster.sort(key = lambda atom: cluster_values[(atom[1], instance.n_values[atom[0]]-1)])

    # output the clusters formatted in LaTex
    res = ''
    for elt in ordered_cluster:
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


#-------------------------------------------------------------------------------
def global_analysis():


    create_regulatory_graph()

    # load the graph from an existing file
    filename = '../../dataset/Imagine/dynamics/regulatory_network_processed.txt'
    graph = Graph('regulatory network')
    graph.load_from_file(filename)

    clusters = {}
    for cluster in graph.clusters:
        if len(graph.clusters[cluster]) >= 20:
            clusters[cluster] = graph.clusters[cluster]
    # print(clusters)

    print('number of clusters in the dynamical graph initially: ', len(graph.clusters))
    print('number of selected clusters (size <= 20): ', len(clusters))

    # load the discrete dataset
    filename = '../../dataset/Imagine/discrete_matrix.csv'
    df_discrete = pd.read_csv(filename, index_col = 0)
    # print(df_discrete.head())

    # load the cell types
    filename = '../../dataset/Imagine/cell_types.csv'
    df_celltypes = pd.read_csv(filename, index_col = 0)
    # print(df_celltypes.head())

    # load the cell macro types
    filename = '../../dataset/Imagine/cell_types_macro.csv'
    df_cell_macrotypes = pd.read_csv(filename, index_col = 0)
    # print(df_cell_macrotypes.head())

    # load the UMAP representation
    filename = '../../dataset/Imagine/umap_coordinates.csv'
    df_umap = pd.read_csv(filename, index_col = 0)
    # print(df_umap.head())

    # load the normalized dataset
    filename = '../../dataset/Imagine/normalized_matrix.csv'
    df_normalized = pd.read_csv(filename, index_col = 0)
    df_normalized = df_normalized.T
    # print(df_normalized.head())

    # fig_comp_graph, axs_comp_graph = plt.subplots(1,2)


    # display the regulatory graph
    print('display the (directed) graph')
    # col_option = '01_colors'
    col_option = 'clustering_colors'
    arrows = True
    cluster_size_limit = 20
    fig, ax = plt.subplots()
    # graph.plot(ax, col_option, arrows, cluster_size_limit)
    ax.set_title('2d representation of the dynamical graph')
    # graph.plot(ax, col_option, arrows, cluster_size_limit, [6,7,8,9,11,20])
    graph.plot(ax, col_option, arrows, cluster_size_limit)
    # ax.set_title('Graphe dynamique')
    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
    ax.set_xlabel('FA 1')
    ax.set_ylabel('FA 2')
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    # col_option = 'clustering_colors'
    # fig, ax = plt.subplots()
    # graph.plot(ax, col_option, arrows, cluster_size_limit)
    # ax.set_title('Dynamical graph')



    # ################################################################################
    # #### load and process the coexpression network
    # # input_file = '../../dataset/Imagine/coexpression/coexpression_network.txt'
    # # n_cell_min = 200
    # # score_min = 0.35
    # # louvain_param = 0.7
    # # output_file = '../../dataset/Imagine/coexpression/coexpression_graph.txt'
    # # codata = GData()
    # # codata.load_from_file(input_file, n_cell_min, score_min)
    # # cograph = Graph('coexpression network')
    # # # build graph from raw data, exclude mitochondrial and ribosomal genes
    # # exclude_mt_rp = True
    # # filter_edges = True
    # # cograph.create_from_gdata(codata, exclude_mt_rp, filter_edges)
    # # cograph.compute_clusters(louvain_param)
    # # cograph.compute_positions()
    # # cograph.save(output_file)
    #
    cograph = Graph('coexpression network')
    cograph.load_from_file('../../dataset/Imagine/coexpression/coexpression_graph.txt')

    # display the coexpression graph
    print('\n\ndisplaying the coexpression graph\n\n')
    col_option = 'clustering_colors'
    fig, ax = plt.subplots()
    # cograph.plot(ax, col_option, False, 10)
    ax.set_title('2d representation of the coexpression graph')
    # cograph.plot(ax, col_option, False, 10, [1,2,5], 'upper right')
    cograph.plot(ax, col_option, False, 10)
    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
    ax.set_xlabel('FA 1')
    ax.set_ylabel('FA 2')
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    print('number of clusters in the coexpression graph initially: ', len(cograph.clusters))
    coexpression_clusters = {}
    for cluster in cograph.clusters:
        if len(cograph.clusters[cluster]) >= 20:
            coexpression_clusters[cluster] = cograph.clusters[cluster]
    # print(coexpression_clusters)

    # fig, ax = plt.subplots()
    # ax.set_title('Coexpression graph')
    # cograph.plot(ax, 'clustering_colors', False, 10)
    #
    # # save the dynamical graph
    # # filename = '../../dataset/Imagine/regulatory_network_processed.txt'
    # # graph.save(filename)





    # ################################################################################
    # #### load the coexpression network on the T cells
    #
    # # analysis of the sub co-expression network on the T cells for comparison
    # filename = '../../dataset/Imagine/coexpression/T_coexpression_network.txt'
    # n_cell_min = 10
    # score_min = 0.35
    # louvain_param = 0.7
    # data_coexpr_T = GData()
    # data_coexpr_T.load_from_file(filename, n_cell_min, score_min)
    # graph_coexpr_T = Graph('coexpression network (T cells)')
    # # build graph from raw data, exclude mitochondrial and ribosomal genes
    # exclude_mt_rp = True
    # filter_edges = True
    # graph_coexpr_T.create_from_gdata(data_coexpr_T, exclude_mt_rp, filter_edges)
    # graph_coexpr_T.compute_clusters(louvain_param)
    # graph_coexpr_T.compute_positions()
    # fig, ax = plt.subplots()
    # ax.set_title('T cells coexpression graph')
    # graph_coexpr_T.plot(ax, 'clustering_colors', False, 10)
    #
    # # plots the clusters
    #
    # clusters_coexpr_T = {}
    # for cluster in graph_coexpr_T.clusters:
    #     if len(graph_coexpr_T.clusters[cluster]) >= 10:
    #         clusters_coexpr_T[cluster] = graph_coexpr_T.clusters[cluster]
    #
    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    #
    # ncol = 2
    # nrows = 2
    # fig, axs = plt.subplots(nrows, ncol)
    # fig.suptitle('clusters graphe coexpr T cells')
    # axs = axs.flat
    # ind_plot = 0
    # for cluster_index in clusters_coexpr_T:
    #     ax = axs[ind_plot]
    #     body =  [ instance.get_atom_index(graph_coexpr_T.atoms[index]) for index in clusters_coexpr_T[ cluster_index ] ]
    #     # fig, ax = plt.subplots()
    #     ax.set_title('cluster ' + str(cluster_index))
    #     plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
    #     ind_plot += 1
    #
    # ##########################################################################






    # ##########################################################################

    print('\n\nplotting the dynamical clusters\n\n')
    # create a fake instance
    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    # cluster_index = list(clusters.keys())[0]

    selected_clusters = [6,7,8,9,11,20]
    # selected_clusters = [cl for cl in clusters]

    ncol = 2
    nrows = int(len(selected_clusters)/float(ncol))
    fig, axs = plt.subplots(nrows, ncol)
    fig.suptitle('Dynamical clusters matching error on the cells')
    axs = axs.flat
    ind_plot = 0
    for cluster_index in selected_clusters:
        ax = axs[ind_plot]
        body =  [ instance.get_atom_index(graph.atoms[index]) for index in clusters[ cluster_index ] ]
        # fig, ax = plt.subplots()
        ax.set_title('cluster ' + str(cluster_index))
        plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
        ind_plot += 1

    # ##########################################################################


    # ##########################################################################

    # print('\n\nplot the coexpression clusters\n\n')
    # # create a fake instance
    # instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    #
    # selected_clusters = [cl for cl in coexpression_clusters]
    #
    # ncol = 2
    # nrows = int(len(selected_clusters)/float(ncol))
    # fig, axs = plt.subplots(nrows, ncol)
    # fig.suptitle('coexpression clusters matching error')
    # axs = axs.flat
    # ind_plot = 0
    # for cluster_index in selected_clusters:
    #     ax = axs[ind_plot]
    #     body =  [ instance.get_atom_index(cograph.atoms[index]) for index in coexpression_clusters[ cluster_index ] ]
    #     # fig, ax = plt.subplots()
    #     ax.set_title('cluster ' + str(cluster_index))
    #     plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
    #     ind_plot += 1

    # ##########################################################################


    # # ##########################################################################
    # # plot cluster 8 and selection of cells
    # fig, ax = plt.subplots()
    # body =  [ instance.get_atom_index(graph.atoms[index]) for index in clusters[ 8 ] ]
    # plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
    # ax.set_title('Cluster 8 error on the cells')
    #
    # # cluster 8 histogram
    # histo = histogram.Histogram(instance, body, histogram.Histogram_type.GLOBAL)
    # values = [ error for error in range(len(histo.positive_histogram)) for _ in histo.positive_histogram[error] ]
    # fig,ax = plt.subplots()
    # ax.hist(values, 50, density=True, edgecolor='black')
    # ax.set_ylabel('proportion de cellules')
    # ax.set_xlabel('erreur de couverture')
    #
    # # threshold <= 140
    # selected_cells = []
    # for i in range(190):
    #     selected_cells += histo.positive_histogram[i]
    #
    # T_cells = list(df_cell_macrotypes.loc[df_cell_macrotypes['cellType_macro'] == 'T'].index)
    #
    # selected_cells = [cell for cell in selected_cells if cell in T_cells]
    # other_cells = [cell for cell in T_cells if not cell in selected_cells]
    #
    # df_selected_cells = pd.DataFrame([1]*len(selected_cells)+[0]*len(other_cells), index = selected_cells+other_cells, columns=['X'])
    # # print(df_selected_cells)
    # # print(selected_cells+other_cells)
    # df_selected_cells.to_csv('selected_cells_dynamics_c20.csv')
    #
    # fig,ax = plt.subplots()
    # ax.scatter(df_umap['UMAP_1'][:], df_umap['UMAP_2'][:], s=5, c=['red'  if cell in selected_cells else 'black' for cell in df_umap.index])
    #
    # # print('cluster 8:')
    # # print_cluster_latex(instance, [graph.atoms[index] for index in clusters[ 8 ]])
    #
    # # plot cluster 20
    # fig, ax = plt.subplots()
    # body =  [ instance.get_atom_index(graph.atoms[index]) for index in clusters[ 20 ] ]
    # plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
    # ax.set_title('Cluster 20 error on the cells')


    # print('cluster 20:')
    # print_cluster_latex(instance, [graph.atoms[index] for index in clusters[ 20 ]])


    # ##########################################################################
    # Output the clusters in csv files
    output_clusters_edges = False
    if output_clusters_edges:

        # create a dataset with all the clusters
        data = []
        for cluster_index in clusters:
            for atom_index in clusters[cluster_index]:
                atom = graph.atoms[atom_index]
                data.append([atom[0], atom[1], cluster_index])

        df_dynamical_clusters = pd.DataFrame(data, columns=['gene', 'gene discrete value', 'cluster index'])
        print(df_dynamical_clusters.head())
        df_dynamical_clusters.to_csv('dynamical_clusters.csv')

        # look for edges connecting the clusters
        connecting_edges = []
        # build a table of the cluster indexes of the atoms
        atom_to_cluster_index = {}
        for cluster_index in clusters:
            for atom_index in clusters[cluster_index]:
                atom_to_cluster_index[graph.atoms[atom_index]] = cluster_index
        for edge in graph.edges:
            atom_left = graph.atoms[edge[0]]
            atom_right = graph.atoms[edge[1]]
            if atom_left in atom_to_cluster_index and atom_right in atom_to_cluster_index:
                if atom_to_cluster_index[atom_left] != atom_to_cluster_index[atom_right]:
                    print('edge: ', atom_left, ' ', atom_to_cluster_index[atom_left], ' -> ', atom_right, ' ', atom_to_cluster_index[atom_right])
                    connecting_edges.append([atom_left[0], atom_left[1], atom_to_cluster_index[atom_left], atom_right[0], atom_right[1], atom_to_cluster_index[atom_right]])

        connecting_edges.sort(key=lambda elt: (elt[2], elt[5], elt[0], elt[3]))

        df_connecting_edges = pd.DataFrame(connecting_edges, columns=['left gene', 'left value', 'left cluster', 'right gene', 'right value', 'right cluster'])
        df_connecting_edges.to_csv('dynmical_connecting_edges.csv')



    ###############################################
    # compare dynamical clusters between themselves
    ###############################################

    paired_clusters = [(0,6), (1,7), (2,8), (3,9), (4,11), (20,23)]

    print('Shared genes between the dynamical clusters:\n')
    for elt in paired_clusters:
        ind_1 = elt[0]
        ind_2 = elt[1]
        genes_c1 = [graph.atoms[index][0] for index in clusters[ind_1]]
        genes_c2 = [graph.atoms[index][0] for index in clusters[ind_2]]
        shared = [elt for elt in genes_c1 if elt in genes_c2]
        print('between cluster ', ind_1, '(', len(genes_c1), ') and ', ind_2, '(', len(genes_c2), ') genes: ', len(shared), ' genes are shared')



    ###############################################
    # compare dynamical and coexpression clusters
    ###############################################

    # compute number of atoms shared between a dynamical and a coexpression cluster
    print('\n\nShared atoms between the dynamical and the coexpression clusters:\n')
    matching_clusters = {}
    for cluster in clusters:
        list_matched = []
        for coexpression_cluster in coexpression_clusters:
            atoms_c1 = [graph.atoms[index] for index in clusters[cluster]]
            atoms_c2 = [cograph.atoms[index] for index in coexpression_clusters[coexpression_cluster]]
            shared = [elt for elt in atoms_c1 if elt in atoms_c2]
            if float(len(shared))/len(atoms_c1) >= 0.1 or float(len(shared))/len(atoms_c2) >= 0.1:
                print('between dynamical cluster ', cluster, '(', len(atoms_c1), ' atoms) and coexpression cluster', coexpression_cluster, '(', len(atoms_c2), ' atoms) : ', len(shared), ' atoms are shared')
                list_matched.append(coexpression_cluster)

        if len(list_matched) > 0:
            matching_clusters[cluster] = list_matched

    # # plot the matched clusters (dynamical and coexpression)
    # for dynamical_index in matching_clusters:
    #
    #     fig, axs = plt.subplots(1, len(matching_clusters[dynamical_index])+1)
    #
    #     body =  [ instance.get_atom_index(graph.atoms[index]) for index in clusters[dynamical_index] ]
    #     plot_cluster_umap(df_discrete, df_umap, instance, body, axs[0])
    #     axs[0].set_title('dynamical cluster ' + str(dynamical_index))
    #
    #     for ind in range(len(matching_clusters[dynamical_index])):
    #         body =  [ instance.get_atom_index(cograph.atoms[index]) for index in coexpression_clusters[matching_clusters[dynamical_index][ind]] ]
    #         plot_cluster_umap(df_discrete, df_umap, instance, body, axs[ind+1])
    #         axs[ind+1].set_title('coexpression cluster ' + str(matching_clusters[dynamical_index][ind]))

    # # plot the atoms that are not associated to coexpression clusters
    # cl = [8,11,20]
    # atoms = []
    # for index in cl:
    #     atoms += [graph.atoms[atom_ind] for atom_ind in clusters[index]]
    # print(atoms)
    # file = open('atoms_unshared.txt', 'w')
    # for atom in atoms:
    #     file.write(atom[0] + ' ' + str(atom[1]) + '\n')
    # file.close()

    # ###############################################
    # # load the old coexpression clusters
    # ###############################################

    # # filename = 'global_clusters_old_1.txt'
    # filename = 'global_clusters_old_2.txt'
    # file = open(filename, 'r')
    # content = file.read().splitlines()
    # file.close()

    # old_clusters_2 = []
    # for line in content:
    #     tokens = line.split(', ')
    #     atoms_ = []
    #     for elts in tokens:
    #         elt = elts.split('_')
    #         atoms_.append((elt[0], int(elt[1])))
    #     old_clusters_2.append(atoms_)
    # print(old_clusters_2)

    # instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    # selected_clusters = old_clusters_2
    # ncol = 4
    # nrows = int(len(selected_clusters)/float(ncol))+1
    # fig, axs = plt.subplots(nrows, ncol)
    # fig.suptitle('coexpression old_1 clusters matching error')
    # axs = axs.flat
    # ind_plot = 0
    # for cluster in selected_clusters:
    #     ax = axs[ind_plot]
    #     body =  [ instance.get_atom_index(atom) for atom in cluster ]
    #     # fig, ax = plt.subplots()
    #     ax.set_title('cluster ' + str(cluster_index))
    #     plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
    #     ind_plot += 1


    # ##########################################################################
    # display the biomarkers from the clusters

    print('\n\nBimarkers in the dynamical clusters:\n\n')

    biom_list = ['AIF1', 'BLK', 'CCL2', 'CCL5', 'CCR6', 'CD14', 'CD16', 'CD180', 'CD19', 'CD37', 'CD3D']
    biom_list += ['CD3E', 'CD40LG', 'CD72', 'CD74', 'CD79A', 'CD79B', 'CD80', 'CD83', 'CD8A', 'CD8B', 'CD96', 'CSF1R', 'CST7']
    biom_list += ['CTLA4', 'CXCR5', 'FASLG', 'FCER1A', 'FCER1G', 'FCGR3A', 'FCN1', 'FCRL2', 'FLT3', 'GNLY', 'GPR183', 'GZMB']
    biom_list += ['GZMH', 'HLA-DMA', 'HLA-DQA1', 'HTR3A', 'IL1R2', 'IL2RB', 'IL32', 'IL6R', 'IL7R', 'ITGAX', 'ITM2A']
    biom_list += ['KIR2DL3', 'KLRC1', 'KLRD1', 'KLRG1', 'LGALS2', 'LGALS3', 'LST1', 'LYZ', 'MS4A1']
    biom_list += ['NCAM1', 'NKG7', 'PAX5', 'PIKFYVE', 'PNOC', 'PRF1', 'S100A4', 'S100A9', 'SPIB', 'STAP1']
    biom_list += ['TBX21', 'TCL1A' , 'TCL1B', 'TLR7', 'TNFRSF13B', 'TYROBP', 'VPREB3', 'ZBTB16']

    for cluster_index in [6,7,8,9,11,20]:
        genes_markers =  [ graph.atoms[index] for index in clusters[ cluster_index ] if graph.atoms[index][0] in biom_list ]

        print('\n\nmarkers for cluster ' + str(cluster_index) + ':')
        # print(genes_markers)

        list_biom_str = ''
        for elt in genes_markers:
            list_biom_str += elt[0] + '_{' + str(elt[1]) + '/' + str(instance.n_values[elt[0]]-1) + '}, '
        print(list_biom_str + '\n')

    # ##########################################################################

    # for ind_cluster in [8, 11, 20]:
    #     cl = clusters[ ind_cluster ].copy()
    #     print('\n\ncluster ', ind_cluster, ' : ', len(cl),  '\n\n\n')
    #     cl.sort(key=lambda index: graph.atoms[index][0])
    #     for index in cl:
    #         print(graph.atoms[index])

    print('\n\nLength of the clusters from the dynamical graph:')
    for ind_cluster in clusters:
        print(len(clusters[ind_cluster]), 'atoms in cluster', ind_cluster, '\n\n')


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

    filename = '../../dataset/Imagine/dynamics/T_regulatory_network.txt'

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
    ax.set_title('T cells regulatory graph')
    graph.plot(ax, 'clustering_colors', True, 10)

    fig, ax = plt.subplots()
    ax.set_title('T cells regulatory graph')
    graph.plot(ax, '01_colors', True, 10)

    # graph.save(output_file)

    plt.show()

    return

#-------------------------------------------------------------------------------
def alternative_network():

    # load the discrete dataset
    filename = '../../dataset/Imagine/discrete_matrix.csv'
    df_discrete = pd.read_csv(filename, index_col = 0)
    print(df_discrete.head())

    # load the cell types
    filename = '../../dataset/Imagine/cell_types.csv'
    df_celltypes = pd.read_csv(filename, index_col = 0)
    print(df_celltypes.head())

    # load the cell macro types
    filename = '../../dataset/Imagine/cell_types_macro.csv'
    df_cell_macrotypes = pd.read_csv(filename, index_col = 0)
    print(df_cell_macrotypes.head())

    # load the UMAP representation
    filename = '../../dataset/Imagine/umap_coordinates.csv'
    df_umap = pd.read_csv(filename, index_col = 0)
    print(df_umap.head())



    ############################################################################
    # load the regulatory graph
    filename = '../../dataset/Imagine/dynamics/dynamical_network_035_me6.txt'

    n_cell_min = 20
    score_min = 0.35
    louvain_param = 0.7
    data = GData()
    data.load_from_file(filename, n_cell_min, score_min)
    graph = Graph('regulatory network (max 6 edges)')
    # build graph from raw data, exclude mitochondrial and ribosomal genes
    exclude_mt_rp = True
    filter_edges = True
    graph.create_from_gdata(data, exclude_mt_rp, filter_edges)
    graph.compute_clusters(louvain_param)
    graph.compute_positions()

    fig, ax = plt.subplots()
    ax.set_title('regulatory graph (max 6 edges)')
    graph.plot(ax, 'clustering_colors', True, 10)

    fig, ax = plt.subplots()
    ax.set_title('regulatory graph (max 6 edges)')
    graph.plot(ax, '01_colors', True, 10)

    output_file = '../../dataset/Imagine/dynamics/dynamical_graph_6me_processed.txt'
    graph.save(output_file)

    clusters = {}
    for cluster in graph.clusters:
        if len(graph.clusters[cluster]) >= 20:
            clusters[cluster] = graph.clusters[cluster]

    print('plot the clusters')

    # create a fake instance
    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    # cluster_index = list(clusters.keys())[0]

    # selected_clusters = [4, 9, 14, 195, 81, 226]
    selected_clusters = [cluster_index for cluster_index in clusters]

    n_clusters =  len(clusters)

    ncol = 3
    nrows = int(math.ceil(len(selected_clusters)/float(ncol)))
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


    ################################################################################


    ############################################################################
    # load the regulatory graph
    filename = '../../dataset/Imagine/coexpression/coexpression_network_me6.txt'

    n_cell_min = 20
    score_min = 0.35
    louvain_param = 0.7
    coexpression_data = GData()
    coexpression_data.load_from_file(filename, n_cell_min, score_min)
    coexpression_graph = Graph('coexpression network (max 6 edges)')
    # build graph from raw data, exclude mitochondrial and ribosomal genes
    exclude_mt_rp = True
    filter_edges = True
    coexpression_graph.create_from_gdata(coexpression_data, exclude_mt_rp, filter_edges)
    coexpression_graph.compute_clusters(louvain_param)
    coexpression_graph.compute_positions()

    fig, ax = plt.subplots()
    ax.set_title('coexpression graph (max 6 edges)')
    coexpression_graph.plot(ax, 'clustering_colors', True, 10)

    fig, ax = plt.subplots()
    ax.set_title('coexpression graph (max 6 edges)')
    coexpression_graph.plot(ax, '01_colors', True, 10)

    # output_file = '../../dataset/Imagine/coexpression/coexpression_graph_6me_processed.txt'
    # coexpression_graph.save(output_file)

    coexpression_clusters = {}
    for cluster in coexpression_graph.clusters:
        if len(coexpression_graph.clusters[cluster]) >= 5:
            coexpression_clusters[cluster] = coexpression_graph.clusters[cluster]

    print('plot the clusters')

    # create a fake instance
    instance = Instance.create_random_instance(df_discrete.copy(deep=False), 0.5)
    # cluster_index = list(clusters.keys())[0]

    # selected_clusters = [4, 9, 14, 195, 81, 226]
    selected_coexpression_clusters = [cluster_index for cluster_index in coexpression_clusters]

    n_clusters =  len(clusters)

    print(len(coexpression_clusters))

    ncol = 3
    nrows = int(math.ceil(len(selected_coexpression_clusters)/float(ncol)))
    fix, axs = plt.subplots(nrows, ncol)
    axs = axs.flat
    ind_plot = 0
    for cluster_index in selected_coexpression_clusters:
        ax = axs[ind_plot]
        body =  [ instance.get_atom_index(coexpression_graph.atoms[index]) for index in coexpression_clusters[ cluster_index ] ]
        # fig, ax = plt.subplots()
        ax.set_title('cluster ' + str(cluster_index))
        plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
        ind_plot += 1

    plt.show()






    return

global_analysis()

# sub_analysis()

# alternative_network()
