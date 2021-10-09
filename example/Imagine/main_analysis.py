#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import discretization

import UMAP_visualisation

import celltype_classification

import network_clustering

import coexpression_analysis


print('************************************************************************************************')
print('    Analysis of a pbmc dataset through quantified inductive logic programming with LOLH')
print('************************************************************************************************')

filename = '../../dataset/normalized_matrix.csv'
filename_discrete = '../../dataset/discrete_dataset.csv'

print('\n\n')
print('1) Discretization of the dataset')

discretization.discretization(filename, filename_discrete)

print('\n\n')
print('2) Visualization of the data (UMAP)')

UMAP_visualisation.UMAP_visualisation()

print('\n\n')
print('3) Cell type classification (NK)')

celltype_classification.celltype_classification()

print('\n\n')
print('4) Computation of a co-expresion network')

print('\n')

cmd = './../../c++/main'
cmd += ' -im ../../dataset/Imagine/discrete_dataset.csv'
cmd += ' -o coexpression_network.txt'
cmd += ' -t 0.3'
os.system(cmd)

print('\n\n')
print('5) Clustering the network')

print('\n')

input_file = './coexpression_network.txt'
n_cell_min = 200
score_min = 0.35
louvain_param = 0.7
output_file = './coexpression/gene_clusters.txt'
network_clustering.network_clustering(input_file, n_cell_min, score_min, louvain_param, output_file)

print('\n\n')
print('6) Analysis of the network')

coexpression_analysis.process_global_clusters()

print('\n\n\n')
print('7) Computation of a sub network')

cmd = './../../c++/main'
cmd += ' -im ../../dataset/Imagine/sub_dataset_discrete.csv'
cmd += ' -o ./sub_coexpression_network.txt'
cmd += ' -t 0.3'
os.system(cmd)

print('\n\n\n')
print('8) Clustering of the sub network')

print('\n')

input_file = './sub_coexpression_network.txt'
n_cell_min = 5
score_min = 0.35
louvain_param = 0.4
output_file = './sub_network_gene_clusters.txt'
network_clustering.network_clustering(input_file, n_cell_min, score_min, louvain_param, output_file)

print('\n\n\n')
print('9) Analysis of the sub network')

coexpression_analysis.process_sub_network_clusters()
