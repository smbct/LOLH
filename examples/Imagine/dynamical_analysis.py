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


# create_regulatory_graph()

# load the graph from an existing file
filename = '../../dataset/Imagine/regulatory_network_processed.txt'
graph = Graph('regulatory network')
graph.load_from_file(filename)


# load the discrete dataset
filename = '../../dataset/Imagine/discrete_matrix.csv'
df = pd.read_csv(filename, index_col = 0)
print(df.head())

# load the cell types
filename = '../../dataset/Imagine/cell_types.csv'
df_celltypes = pd.read_csv(filename, index_col = 0)
print(df_celltypes.head())

# load the UMAP representation
filename = '../../dataset/Imagine/umap_coordinates.csv'
df_umap = pd.read_csv(filename, index_col = 0)
print(df_umap.head())

# load the normalized dataset

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

# save the graph
# filename = '../../dataset/Imagine/regulatory_network_processed.txt'
# graph.save(filename)

# clusters = {}
# for cluster in graph.clusters:
#     if len(graph.clusters[cluster]) >= 20:
#         clusters[cluster] = graph.clusters[cluster]
# print(clusters)

plt.show()
