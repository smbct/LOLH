#!/usr/bin/python

import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import os

import sys
sys.path.append('../../python')

from instance import Instance

import visualizer
import histogram

from network import GData
from network import Graph

# # read the atoms from the file
# filename = 'atoms_unshared.txt'
# file = open(filename, 'r')
# content = file.read().splitlines()
# file.close()
# atoms = []
# for elt in content:
#     tokens = elt.split(' ')
#     atoms.append((tokens[0], int(tokens[1])))
#
# # read the coexpression graph
# filename = '../../dataset/Imagine/coexpression/coexpression_network.txt'
# file = open(filename, 'r')
# content = file.read().splitlines()
# file.close()
#
# content_unshared = []
#
# for line in content:
#     tokens = line.split(' ')
#     atom = (tokens[0], int(tokens[1]))
#     if atom in atoms:
#         content_unshared.append(line)
# print(len(content_unshared))
#
# # save the sub network
# filename = '../../dataset/Imagine/coexpression/coexpression_network_sub_test_dynamics.txt'
# file = open(filename, 'w')
# ind = 0
# for line in content_unshared:
#     file.write(line)
#     if ind < len(content_unshared) - 1:
#         file.write('\n')
#         ind += 1
# file.close()


# create a sub coexpression graph

input_file = '../../dataset/Imagine/coexpression/coexpression_network_sub_test_dynamics.txt'
n_cell_min = 1
score_min = 0.005
louvain_param = 0.7
# output_file = '../../dataset/Imagine/coexpression/coexpression_graph_sub_test_dynamics.txt'
data = GData()
data.load_from_file(input_file, n_cell_min, score_min)
graph = Graph('sub coexpression network')
# build graph from raw data, exclude mitochondrial and ribosomal genes
exclude_mt_rp = True
filter_edges = True
graph.create_from_gdata(data, exclude_mt_rp, filter_edges)
graph.compute_clusters(louvain_param)
graph.compute_positions()
# graph.save(output_file)

# display the coexpression graph
print('display the coexpression graph')
col_option = 'clustering_colors'
fig, ax = plt.subplots()
graph.plot(ax, col_option, False)
ax.set_title('Coexpression graph')

plt.show()
