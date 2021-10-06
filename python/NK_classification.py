#!/usr/bin/python3

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.patches import Patch

# LOLH objects
from instance import Instance
from solver import Solver

import visualizer
import histogram

# import the cell types
file_name = '../dataset/Imagine/cell_types.csv'
df_cell_types = pd.read_csv(file_name, index_col=0)
df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)

# import the UMAP 2d representation
file_name = '../dataset/Imagine/embedding_coord.csv'
umap_coordinates = pd.read_csv(file_name, index_col = 0)

# import the single cell discretized matrix
file_name = '../dataset/Imagine/IMAGINE_normalised_discrete_adaptive.csv'
df = pd.read_csv(file_name, index_col=0)

# creation of the classification instance
cell_type = 'NK'
instance = Instance.create_cluster_instance(df, df_cell_types, cell_type)

positive_cells = instance._pos_samples
negative_cells = instance._neg_samples

# print('n atoms: ', inst.n_atoms())
print('n positive samples: ', instance.n_positives())
print('n negative samples: ', instance.n_negatives())

# display the atom scores

# selected_atoms, scores = solver.select_best_atoms_threshold(0.55)

fig, ax = plt.subplots()
ax.scatter([instance.atom_score[atom_index][0] for atom_index in instance.n_atoms()], [instance.atom_score[atom_index][1] for atom_index in instance.n_atoms()], marker='x')
ax.plot([0, instance.n_positives()], [0, instance.n_negatives()])

plt.show()
