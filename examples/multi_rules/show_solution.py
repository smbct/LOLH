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


print('\n')
print('Classification of the NK cells')
print('\n')

# read the discrete matrix
filename = '../../dataset/Imagine/discrete_matrix.csv'
df_discrete = pd.read_csv(filename, index_col=0)
# print(df_discrete.head())


# read the normalized matrix
# filename = '../../dataset/Imagine/normalized_matrix.csv'
# df_normalized = pd.read_csv(filename, index_col=0).T
# print(df_normalized.head())

# read the UMAP coordinates
filename = '../../dataset/Imagine/umap_coordinates.csv'
df_umap = pd.read_csv(filename, index_col = 0)

# read the seurat cell types
df_celltypes = pd.read_csv('../../dataset/Imagine/cell_types.csv', index_col = 0)
df_celltypes.rename(columns={'cellType_final': 'Label'}, inplace=True)
# print(df_celltypes.head())

df_macrotypes = pd.read_csv('../../dataset/Imagine/cell_types_macro.csv', index_col = 0)
df_macrotypes.rename(columns={'cellType_macro': 'Label'}, inplace=True)

# initialization of the classification instance: classification of the NK cells
celltype = 'NK'
celltype = 'T'
instance = Instance.create_cluster_instance(df_discrete.copy(deep=False), df_macrotypes, celltype)

print('Classification of the NK cells')
print('- ', instance.n_positives(), ' positive examples')
print('- ', instance.n_negatives(), ' negative examples')
print('\n')

# load the solution
file = open('T_sol.txt', 'r')
content = file.read().split(' ')
content = content[1:]
sol_values = [int(elt) for elt in content]
print(sol_values)

solver = Solver(instance)

# create a rule from a threshold
sorted_atoms, scores = solver.select_best_atoms_fast()
threshold = 0.5

pred_col = ['blue', 'green', 'orange']
colors = {}
for elt in df_umap.index:
    colors[elt] = 'red'
for ind in range(len(sol_values)):
    val = sol_values[ind]
    colors[instance._pos_samples[ind]] = pred_col[val]

fig, ax = plt.subplots()
ax.scatter(df_umap['UMAP_1'][:], df_umap['UMAP_2'][:], marker='o', s=1, c=[colors[index] for index in df_umap.index])
ax.set_aspect((ax.get_xlim()[1] - ax.get_xlim()[0])/(ax.get_ylim()[1] - ax.get_ylim()[0]))

plt.show()
