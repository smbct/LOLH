#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import random

import os

import umap_visualisation

import NK_classification

import network_clustering

import coexpression_analysis


# fix the random seed
np.random.seed(42)
random.seed(42)

print('*******************************************************************************')
print('  Coexpression analysis of the pbmc dataset from Imagine institute with LOLH   ')
print('*******************************************************************************')

filename = '../../dataset/Imagine/normalized_matrix.csv'
filename_discrete = '../../dataset/discrete_matrix.csv'


print('\n\n')
print('1) Visualization of the data (UMAP)')

umap_visualisation.umap_visualisation()

print('\n\n')
print('2) Cell type classification (NK)')

NK_classification.NK_classification()

print('\n\n')
print('3) Computation of a co-expresion network')

print('\n')

cmd = './../../c++/main'
cmd += ' -im ../../dataset/Imagine/discrete_matrix.csv'
cmd += ' -o ../../dataset/Imagine/coexpression/coexpression_network.txt'
cmd += ' -t 0.3'
os.system(cmd)

print('\n\n')
print('4) Clustering the network')

print('\n')

input_file = '../../dataset/Imagine/coexpression/coexpression_network.txt'
n_cell_min = 200
score_min = 0.35
louvain_param = 0.7
output_file = '../../dataset/Imagine/coexpression/gene_clusters.txt'
network_clustering.network_clustering(input_file, n_cell_min, score_min, louvain_param, output_file)

print('\n\n')
print('5) Analysis of the network')

coexpression_analysis.process_global_clusters()

print('\n\n\n')
print('6) Computation of a sub network')

cmd = './../../c++/main'
cmd += ' -im ../../dataset/Imagine/sub_dataset_discrete.csv'
cmd += ' -o ../../dataset/Imagine/coexpression/sub_coexpression_network.txt'
cmd += ' -t 0.3'
os.system(cmd)

print('\n\n\n')
print('7) Clustering of the sub network')

print('\n')

input_file = '../../dataset/Imagine/coexpression/sub_coexpression_network.txt'
n_cell_min = 5
score_min = 0.35
louvain_param = 0.7

output_file = '../../dataset/Imagine/coexpression/sub_network_gene_clusters.txt'
network_clustering.network_clustering(input_file, n_cell_min, score_min, louvain_param, output_file)

print('\n\n\n')
print('8) Analysis of the sub network')

coexpression_analysis.process_sub_network_clusters()
