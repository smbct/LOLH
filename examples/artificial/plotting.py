#!/usr/bin/python

import numpy as np
import pandas as pd
import random

# dimensionality reduction in the data
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib import cm


# LOLH objects
import sys
sys.path.append('../../python')

from instance import Instance
from solver import Solver
import histogram
import visualizer



#------------------------------------------------------------------------------
def compute_A_B(instance, score):
    return (0, score*instance.n_negatives()), ((1.-score)*instance.n_positives(), instance.n_negatives())



random_seed = 42
random.seed(random_seed)
np.random.seed(random_seed)


# load the matrix
df = pd.read_csv('../../dataset/artificial/artificial_matrix.csv', index_col=0)



# classification from positive vs negative
pos_samples = ['s_'+str(ind) for ind in range(0, int(df.shape[0]/2.))]
neg_samples = ['s_'+str(ind) for ind in range(int(df.shape[0]/2.)+1, df.shape[0])]

instance_name = 'artificial'

instance = Instance.create_instance_explicit(df, pos_samples, neg_samples)

# extract positive/negative solutions for solution visualisation
positive_cells = instance._pos_samples
negative_cells = instance._neg_samples

print('Original instance:')
# print('n atoms: ', inst.n_atoms())
print('n positive samples: ', instance.n_positives())
print('n negative samples: ', instance.n_negatives())


# LOLH solver
solver = Solver(instance)




# PCA projection
fig,ax = plt.subplots()
X = df.values.copy()
X=StandardScaler().fit_transform(X)
X_pca = PCA(n_components=2).fit_transform(X)
col = ['forestgreen' for _ in range(int(df.shape[0]/2))] + ['darkred' for _ in range(int(df.shape[0]/2))]
ax.scatter(X_pca[:,0], X_pca[:,1], c=col, s=1)
# ax.set_title('projection ACP')
ax.set_xlabel('ACP 1')
ax.set_ylabel('ACP 2')
ax.legend(loc='upper right', handles=[mpatches.Patch(color='forestgreen', label='exemples positifs'), mpatches.Patch(color='darkred', label='exemples négatifs')])
ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))


# plot atoms errors
fig, ax = plt.subplots()
ax.set_xlim((0,instance.n_positives()))
ax.set_ylim((0,instance.n_negatives()))
ax.set_xlabel('erreur positive')
ax.set_ylabel('erreur négative')


col = [instance.atom_score[ind_atom][1]/instance.n_negatives() - instance.atom_score[ind_atom][0]/instance.n_positives() for ind_atom in range(instance.n_atoms())]
ax.scatter([instance.atom_score[ind_atom][0] for ind_atom in range(instance.n_atoms())], [instance.atom_score[ind_atom][1] for ind_atom in range(instance.n_atoms())], alpha=1, marker='x', zorder=0, label='atomes logiques', c=col)

ax.plot([0, instance.n_positives()], [0, instance.n_negatives()], color='firebrick', label='score = 0.0', lw=1.5, zorder=1)

A,B = compute_A_B(instance, 0.5)
ax.plot([A[0], B[0]], [A[1], B[1]], color='darkviolet', linestyle='--', label='score = 0.5', lw=1.5, zorder=1)

A,B = compute_A_B(instance, -0.5)
ax.plot([A[0], B[0]], [A[1], B[1]], color='goldenrod', linestyle=(0, (5, 1)), label='score = -0.5', lw=1.5, zorder=1)

viridis = cm.get_cmap('viridis', 256)
cnorm = mcolors.Normalize(vmin=-1, vmax=1)
cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=cnorm, cmap=viridis), ax=ax)
cbar.set_label('score')

ax.legend(loc='lower left', bbox_to_anchor= (0.01, 0.01))
ax.set_aspect('equal')
fig.tight_layout()

plt.show()
