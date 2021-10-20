#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/python3

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from matplotlib import cm
from matplotlib import colors as mcolors

# LOLH objects

import sys
sys.path.append('../../python')

from instance import Instance
from solver import Solver

import visualizer
import histogram

#-------------------------------------------------------------------------------
def compute_polygons(instance, nondominated_atoms):

    # pre: ordered non dominated atoms (indexes)

    points_inside = []
    points_outside = []

    index = 0
    for atom_index in nondominated_atoms:
        score = instance.atom_score[atom_index]

        if index == 0:
            bottom = 0
        else:
            bottom = instance.atom_score[nondominated_atoms[index-1]][1]

        if index == len(nondominated_atoms)-1:
            right = instance.n_positives()
        else:
            right = instance.atom_score[nondominated_atoms[index+1]][0]

        points_inside.append((score[0], bottom))
        points_inside.append(score)
        points_inside.append((right, score[1]))

        points_outside.append((score[0], bottom))
        points_outside.append(score)
        points_outside.append((right, score[1]))

        index += 1

    points_inside.append((instance.n_positives(), 0))
    points_inside.append((0,0))

    points_outside.append((instance.n_positives(), instance.n_negatives()))
    points_outside.append((0, instance.n_negatives()))

    return points_inside, points_outside



#-------------------------------------------------------------------------------
def plot_nondominated_atoms(ax, instance):

    # compute the non dominated atoms
    solver = Solver(instance)
    nondominated_atoms = solver.compute_nondominated_atoms(list(range(instance.n_atoms())))
    print('non dominated atoms: ', nondominated_atoms)

    other_atoms_indexes = [atom_index for atom_index in range(instance.n_atoms()) if not atom_index in nondominated_atoms]

    # ax.scatter([instance.atom_score[atom_index][0] for atom_index in LOLH_rule], [instance.atom_score[atom_index][1] for atom_index in LOLH_rule], marker='o', color='red', s=4, label='atomes sélectionnés', zorder=2)

    ax.scatter([instance.atom_score[atom_index][0] for atom_index in nondominated_atoms], [instance.atom_score[atom_index][1] for atom_index in nondominated_atoms], marker='x', color='red', label='atomes non dominés', zorder=2)

    ax.scatter([instance.atom_score[atom_index][0] for atom_index in other_atoms_indexes], [instance.atom_score[atom_index][1] for atom_index in other_atoms_indexes], marker='x', zorder=1, label='atomes dominés', alpha=0.7)

    points_inside, points_outside = compute_polygons(instance, nondominated_atoms)

    # display the dominance cones
    index = 0
    for atom_index in nondominated_atoms:
        score = instance.atom_score[atom_index]
        if index == 0:
            bottom = 0
        else:
            bottom = instance.atom_score[nondominated_atoms[index-1]][1]
        if index == len(nondominated_atoms)-1:
            right = instance.n_positives()
        else:
            right = instance.atom_score[nondominated_atoms[index+1]][0]

        # bottom
        ax.plot([score[0], score[0]], [score[1], bottom], '--', zorder=1, color='black')

        # right
        ax.plot([score[0], right], [score[1], score[1]], '--', zorder=1, color='black')

        # label
        # ax.text(score[0], score[1], str(ind))
        index += 1

    # equation score = 0
    # ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='black', label='score=0.0', zorder=3)

    ax.set_xlim([0, instance.n_positives()])
    ax.set_ylim([0, instance.n_negatives()])

    ax.set_xlabel('erreur positive')
    ax.set_ylabel('erreur négative')
    # ax.set_title('Atom errors for the classification of ' + str(cell_type))

    ax.set_aspect(instance.n_positives()/instance.n_negatives())

    ax.fill([point[0] for point in points_outside], [point[1] for point in points_outside], alpha=0.3, zorder=0, color='red')
    ax.fill([point[0] for point in points_inside], [point[1] for point in points_inside], alpha=0.3, zorder=0, color='green')

    ax.legend(loc='lower right')

    return

# import the cell types
file_name = '../../dataset/Imagine/cell_types.csv'
df_cell_types = pd.read_csv(file_name, index_col=0)
df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)

# import the UMAP 2d representation
file_name = '../../dataset/Imagine/umap_coordinates.csv'
umap_coordinates = pd.read_csv(file_name, index_col = 0)

# import the single cell discretized matrix
file_name = '../../dataset/Imagine/discrete_matrix.csv'
df = pd.read_csv(file_name, index_col=0)


# Measure of the quality for the CD8 instance

# creation of the classification instance for CD8
instance_CD8 = Instance.create_cluster_instance(df.copy(deep=False), df_cell_types, 'CD8')
fig, ax = plt.subplots()
plot_nondominated_atoms(ax, instance_CD8)
solver_CD8 = Solver(instance_CD8)
print('relative area for CD8: ', solver_CD8.compute_relative_atom_area(solver_CD8.compute_nondominated_atoms(list(range(instance_CD8.n_atoms())))))

# Measure of the quality for the NK and CD4 instance
instance_NK = Instance.create_cluster_instance(df.copy(deep=False), df_cell_types, 'NK')
instance_CD4 = Instance.create_cluster_instance(df.copy(deep=False), df_cell_types, 'CD4')

solver_NK = Solver(instance_NK)
solver_CD4 = Solver(instance_CD4)
print('relative area for NK: ', solver_NK.compute_relative_atom_area(solver_NK.compute_nondominated_atoms(list(range(instance_NK.n_atoms())))))
print('relative area for CD4: ', solver_CD4.compute_relative_atom_area(solver_CD4.compute_nondominated_atoms(list(range(instance_CD4.n_atoms())))))

fig, axs = plt.subplots(1, 2)
plot_nondominated_atoms(axs[0], instance_NK)
plot_nondominated_atoms(axs[1], instance_CD4)
axs[0].set_title('Instance NK')
axs[1].set_title('Instance CD4')

plt.show()
