#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import random

import os

import transitions_parameters_analysis

import dynamical_analysis

# fix the random seed
np.random.seed(42)
random.seed(42)

print('****************************************************************************')
print('  Dynamical analysis of the pbmc dataset from Imagine institute with LOLH   ')
print('****************************************************************************')


print('\n\n')
print('1) Extraction of transitions from the neighborhood graph with Seurat')

os.system('Rscript extract_transitions.R')

print('\n\n')
print('2) Computation of atom correlations for a set of parameters')

cmd = './../../c++/main'
cmd += ' -r -cq'
cmd += ' -im ../../dataset/Imagine/discrete_matrix.csv'
cmd += ' -o transitions_quality.txt'
cmd += ' -it ../../dataset/Imagine/transitions.csv'
cmd += ' -tr 0.6 -pnq 0 -td 1 -v'
os.system(cmd)

print('\n\n')
print('3) Analysis of the correlations for a set of parameters')

transitions_parameters_analysis.selected_parameters_quality()

print('\n\n')
print('3) Computation of a dynamical interaction network')

cmd = './../../c++/main'
cmd += ' -r'
cmd += ' -im ../../dataset/Imagine/discrete_matrix.csv'
cmd += ' -o ../../dataset/Imagine/dynamics/regulatory_network03.txt'
cmd += ' -t 0.3'
cmd += ' -it ../../dataset/Imagine/transitions.csv'
cmd += ' -tr 0.6'
cmd += ' -pnq 0 -td 1 -v'
os.system(cmd)


print('\n\n')
print('4) Analysis of the dynamical network')

dynamical_analysis.global_analysis()



###################
# ./main -r -cq -im "../dataset/Imagine/discrete_matrix.csv" -o "transitions_quality.txt" -t 0.3 -it "../dataset/Imagine/transitions.csv" -tr 0.5 -pnq 0 -td 1 -v


# ./main -cq -im "../dataset/Imagine/discrete_matrix.csv" -o "coexpression_quality.txt" -v
