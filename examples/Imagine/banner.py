#!/usr/bin/python

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors as mcolors
import matplotlib.ticker as ticker


import sys
sys.path.append('../../python')

from instance import Instance
from solver import Solver
import visualizer
import histogram

from network import Graph


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
    ax.scatter(df_coordinates_sorted['UMAP_1'].values, df_coordinates_sorted['UMAP_2'].values, c=col, cmap=plasma, norm=cnorm, s=2)

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

# read the discrete matrix
filename = '../../dataset/Imagine/discrete_matrix.csv'
df_discrete = pd.read_csv(filename, index_col=0)
# print(df_discrete.head())

# read the normalized matrix
filename = '../../dataset/Imagine/normalized_matrix.csv'
df_normalized = pd.read_csv(filename, index_col=0).T
# print(df_normalized.head())

# read the seurat cell types
df_celltypes = pd.read_csv('../../dataset/Imagine/cell_types.csv', index_col = 0)
df_celltypes.rename(columns={'cellType_final': 'Label'}, inplace=True)
# print(df_celltypes.head())

# read the 2d cell types
filename = '../../dataset/Imagine/umap_coordinates.csv'
df_umap = pd.read_csv(filename, index_col = 0)

# load the graph from an existing file
filename = '../../dataset/Imagine/dynamics/regulatory_network_processed.txt'
graph = Graph('regulatory network')
graph.load_from_file(filename)

fig, axs = plt.subplots(2,3)

axs[1][0].remove()
axs[1][2].remove()


# initialization of the classification instance: classification of the NK cells
celltype = 'CD8'
instance = Instance.create_cluster_instance(df_discrete.copy(deep=False), df_celltypes, celltype)

print('Classification of the NK cells')
print('- ', instance.n_positives(), ' positive examples')
print('- ', instance.n_negatives(), ' negative examples')
print('\n')

solver = Solver(instance)



# create a rule from a threshold
sorted_atoms, scores = solver.select_best_atoms_fast()
threshold = 0.45
# print('threshold: ', threshold)
# print(scores[:10])
threshold *= instance.n_positives()*instance.n_negatives()
body = [sorted_atoms[ind] for ind in range(len(sorted_atoms)) if scores[ind] >= threshold] # 21 atoms selected

print('-> Output of LOLH: ', [instance.get_atom(index)[0] + '_' + str(instance.get_atom(index)[1]) for index in body])



# visualize the macthing error on the UMAP
ax = axs[0][0]
plot_cluster_umap(df_discrete, df_umap, instance, body, ax)
ax.set_title('CD8 instance matching error')



# create an histogram
histo = histogram.Histogram(instance, body)
ax = axs[0][1]
visualizer.plot_histograms(ax, histo)
ax.set_title('LOLH rule matching error for the NK cells')
ax.set_ylim([0, 0.3])




# scatter plot of the atoms for the NK instance
celltype = 'NK'
instance = Instance.create_cluster_instance(df_discrete.copy(deep=False), df_celltypes, celltype)
sorted_atoms, scores = solver.select_best_atoms_fast()
threshold = 0.55
# body = [sorted_atoms[ind] for ind in range(len(sorted_atoms)) if scores[ind] >= threshold*instance.n_positives()*instance.n_negatives()] # 21 atoms selected
body = [ind for ind in range(instance.n_atoms()) if instance.atom_score[ind][1]/instance.n_negatives() - instance.atom_score[ind][0]/instance.n_positives() >= threshold]
print('NK instance solution: ')
print([instance.get_atom(ind) for ind in body])
other_atoms = [ind for ind in range(instance.n_atoms()) if not ind in body]
ax = axs[0][2]
ax.set_xlim((0, instance.n_positives()))
ax.set_ylim((0, instance.n_negatives()))

for ind_atom in range(instance.n_atoms()):
    pos = instance.atom_score[ind_atom]
    score = pos[1]/instance.n_negatives() - pos[0]/instance.n_positives()
    atom = instance.get_atom(ind_atom)
    if score >= 0.6:
        ax.text(pos[0], pos[1], '$' + atom[0] + '_{' + str(atom[1]) + '/' + str(instance.n_values[atom[0]]-1) + '}$', zorder = 3).set_clip_on(True)

ax.scatter([instance.atom_score[index][0] for index in body], [instance.atom_score[index][1] for index in body], marker='x', label='selected atoms', c='red', zorder = 1)
ax.scatter([instance.atom_score[index][0] for index in other_atoms], [instance.atom_score[index][1] for index in other_atoms], marker='x', label='other atoms', c='C0', zorder = 1)

ax.plot([0,instance.n_positives()], [0,instance.n_negatives()], '--', c='k', label='score = 0.0', zorder=2)
ax.plot([0, (1-threshold)*instance.n_positives()], [threshold*instance.n_negatives(), instance.n_negatives()], '--', c='red', label='score = '+str(threshold), zorder=2)

ax.set_aspect( (ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]) )
ax.set_xlabel('positive error')
ax.set_ylabel('negative error')
ax.set_title('atom errors for NK classification')
ax.legend(loc='lower left')




# plot the dynamical graph
col_option = 'clustering_colors'
arrows = True
cluster_size_limit = 20
ax = axs[1][1]
ax.set_title('2d representation of the dynamical graph')
graph.plot(ax, col_option, arrows, cluster_size_limit, [6,7,8,9,11,20])
# ax.set_title('Graphe dynamique')
# ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
# remove FA 2 coordinates
ax.xaxis.set_major_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.NullLocator())
ax.set_xlabel('FA 1')
ax.set_ylabel('FA 2')




plt.show()
