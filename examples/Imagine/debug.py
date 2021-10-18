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

# creation of the classification instance
cell_type = 'NK'
instance = Instance.create_cluster_instance(df.copy(deep=False), df_cell_types, cell_type)

positive_cells = instance._pos_samples
negative_cells = instance._neg_samples

score_thresold = 0.55

# compute the body from the score threshold
solver = Solver(instance)
LOLH_rule, _ = solver.select_best_atoms_threshold(score_thresold)

print([instance.get_atom(index) for index in LOLH_rule])
