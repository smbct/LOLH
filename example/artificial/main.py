#!/usr/bin/python

import os

import generator
import evaluation

print('********************************************')
print('****** Generation of a random dataset ******')
print('********************************************')

os.system('mkdir ../../dataset/artificial')

filename = '../../dataset/artificial/artificial_matrix.csv'
nVariables = 10000
nPositives = 1000
nNegatives = 1000

generator.generate_dataset(nVariables, nPositives, nNegatives, filename)
