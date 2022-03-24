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
# celltype = 'NK'
# instance = Instance.create_cluster_instance(df_discrete.copy(deep=False), df_celltypes, celltype)

celltype = 'T'
instance = Instance.create_cluster_instance(df_discrete.copy(deep=False), df_macrotypes, celltype)


# export the instance to a file
file = open('T_instance.txt', 'w')
file.write(str(len(instance._pos_samples)) + ' ')
for ind in range(len(instance._pos_samples)):
    file.write(instance._pos_samples[ind])
    if ind < len(instance._pos_samples)-1:
        file.write(' ')
file.write('\n')

file.write(str(len(instance._neg_samples)) + ' ')
for ind in range(len(instance._neg_samples)):
    file.write(instance._neg_samples[ind])
    if ind < len(instance._neg_samples)-1:
        file.write(' ')
file.close()

print('Classification of the NK cells')
print('- ', instance.n_positives(), ' positive examples')
print('- ', instance.n_negatives(), ' negative examples')
print('\n')

solver = Solver(instance)

# create a rule from a threshold
sorted_atoms, scores = solver.select_best_atoms_fast()
threshold = 0.55

threshold *= instance.n_positives()*instance.n_negatives()
lolh_body = [sorted_atoms[ind] for ind in range(len(sorted_atoms)) if scores[ind] >= threshold] # 21 atoms selected

print('-> Output of LOLH: ', [instance.get_atom(index)[0] + '_' + str(instance.get_atom(index)[1]) for index in lolh_body])

print('Visualisations of the LOLH rule')

print('\n')

# create an histogram
histo = histogram.Histogram(instance, lolh_body)
fig, ax = plt.subplots()
visualizer.plot_histograms(ax, histo)
ax.set_title('LOLH rule matching error for the NK cells')
ax.set_ylim([0, 0.3])

# plot the UMAP with the matching error
# fig, ax = plt.subplots()

plasma = cm.get_cmap('plasma_r', 256)
cnorm = mcolors.Normalize(vmin=0, vmax=len(histo.positive_histogram))

histo = histogram.Histogram(instance, lolh_body, histogram.Histogram_type.GLOBAL)
cell_score = {barcode : error for error in range(len(histo.positive_histogram)) for barcode in histo.positive_histogram[error]}


# sort all the cells according to the score
barcodes = [index for index in df_umap.index]
barcodes.sort(key = lambda elt: cell_score[elt], reverse=True)
df_coordinates_sorted = df_umap.loc[barcodes]

col = [cell_score[barcode] for barcode in df_coordinates_sorted.index]
fig, ax = plt.subplots()
ax.scatter(df_coordinates_sorted['UMAP_1'][:], df_coordinates_sorted['UMAP_2'][:], marker='o', s=1, c=col, cmap=plasma, norm=cnorm)
ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')
cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=plasma), ax=ax)
cbar.set_label('Matching error')

plt.show()
