#!/usr/bin/python

import numpy as np
import math

from instance import Instance

from enum import Enum

#----------------------------------------------------------------------
class Histogram_type(Enum):
    # numerating value to define histogram type
    POSITIVE = 0
    NEGATIVE = 1
    POSITIVE_NEGATIVE = 2
    GLOBAL = 3

# Define a Histogram class to represent the error histogram of a rule body over the dataset
#----------------------------------------------------------------------
class Histogram:


    #---------------------------------------------------------------------------
    def __init__(self, instance, body, type=Histogram_type.POSITIVE_NEGATIVE):

        self.instance = instance
        self.body = body
        self.type = type

        self.positive_histogram = None
        self.negative_histogram = None

        self._build_histogram()

        return


    #---------------------------------------------------------------------------
    def _build_histogram(self):

        # print('histogram built here')

        body_features = [self.instance.atoms[index][0] for index in self.body]
        
        # in case there is a duplicated gene in the atoms: ex A_0 and A_1 (otherwise there would be duplicated columns by panda)
        body_features = np.unique(body_features)

        if self.type != Histogram_type.GLOBAL:
            if self.type != Histogram_type.NEGATIVE:
                self._build_positive(body_features)
    
            if self.type != Histogram_type.POSITIVE:
                self._build_negative(body_features)
        else:
            self._build_global(body_features)
            
        # print([len(elt) for elt in self.positive_histogram])
        # print([len(elt) for elt in self.negative_histogram])

        return

    #---------------------------------------------------------------------------
    def _build_positive(self, body_features):
        self.positive_histogram = [[] for _ in range(len(self.body)+1)]
        reduced_pos_dataset = self.instance.dataset.loc[self.instance._pos_samples,body_features]

        for pos_index in self.instance._pos_samples:
            sample = reduced_pos_dataset.loc[pos_index]
            error = self._rule_error(sample)
            self.positive_histogram[error].append(pos_index)
        return

    #---------------------------------------------------------------------------
    def _build_negative(self, body_features):
        self.negative_histogram = [[] for _ in range(len(self.body)+1)]
        reduced_neg_dataset = self.instance.dataset.loc[self.instance._neg_samples,body_features]
        for neg_index in self.instance._neg_samples:
            sample = reduced_neg_dataset.loc[neg_index]
            error = self._rule_error(sample)
            self.negative_histogram[error].append(neg_index)
        return
    
    #---------------------------------------------------------------------------
    def _build_global(self, body_features):
        # global is stored into positive histogram only
        self.positive_histogram = [[] for _ in range(len(self.body)+1)]
        reduced_dataset = self.instance.dataset.loc[:,body_features]
        for index in self.instance.dataset.index:
            sample = reduced_dataset.loc[index]
            error = self._rule_error(sample)
            self.positive_histogram[error].append(index)
        return

    #---------------------------------------------------------------------------
    def _rule_error(self, sample) -> int:
        error = 0
        for atom_index in self.body:
            atom = self.instance.atoms[atom_index]
            if sample[atom[0]] != atom[1]:
                error += 1
        return error

    #---------------------------------------------------------------------------
    def positive_covered(self, threshold) -> int:

        if self.type == Histogram_type.NEGATIVE:
            raise Exception('Positive samples not supported in this histogram')

        if threshold >= len(self.positive_histogram):
            raise Exception('Threshold is too high')

        # returns the number of samples covered from the positive set: error <= threshold
        rate = np.sum([len(elt) for elt in self.positive_histogram[:threshold+1]])
        return rate

    #---------------------------------------------------------------------------
    def negative_covered(self, threshold) -> int:

        if self.type == Histogram_type.POSITIVE:
            raise Exception('Negative samples not supported in this histogram')

        if threshold >= len(self.negative_histogram):
            raise Exception('Threshold is too high')

        # returns the number of samples covered from the positive set: error <= threshold
        rate = np.sum([len(elt) for elt in self.negative_histogram[:threshold+1]])
        return rate
