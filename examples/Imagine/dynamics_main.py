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

print('****************************************************************************')
print('  Dynamical analysis of the pbmc dataset from Imagine institute with LOLH   ')
print('****************************************************************************')

filename = '../../dataset/Imagine/normalized_matrix.csv'
filename_discrete = '../../dataset/discrete_matrix.csv'

print('\n\n')
print('1) Computation of a dynamical interaction network')

# run R script automatically
# Rscript extract_transitions.R

# ./main -r -cq -im "../dataset/Imagine/discrete_matrix.csv" -o "transitions_quality.txt" -t 0.3 -it "../dataset/Imagine/transitions.csv" -tr 0.5 -pnq 0 -td 1 -v


# ./main -cq -im "../dataset/Imagine/discrete_matrix.csv" -o "coexpression_quality.txt" -v
