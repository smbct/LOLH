#!/usr/bin/python

import numpy as np
import math

from instance import Instance

# Define a Rule class to represent a logic rule
#----------------------------------------------------------------------
class Rule:

    #---------------------------------------------------------------------------
    def __init__(self, instance, body = []):

        self.atomList = atomList
        self.body = body
        self.score = (0,0)
        self.positive_histogram = []
        self.negative_histogram = []

        self.pos_error_indexes = [[] for _ in range(len(body)+1)]
        self.neg_error_indexes = [[] for _ in range(len(body)+1)]

        if len(self.body) > 0:
            self.init_body(body)

        return


    #---------------------------------------------------------------------------
    def init_body(self, body):
        self.body = body.copy()
        self.compute_score()
        self.compute_histograms()
        return

    #---------------------------------------------------------------------------
    def compute_score(self):
        pos_sum = 0
        neg_sum = 0
        for ind in self.body:
            pos_sum += self.atomList.score[ind,0]
            neg_sum += self.atomList.score[ind,1]
        self.score = (pos_sum, neg_sum)
        return

    #---------------------------------------------------------------------------
    def compute_histograms(self):
        self.positive_histogram = [0 for _ in range(len(self.body)+1)]
        self.negative_histogram = [0 for _ in range(len(self.body)+1)]

        # matching score on positive samples
        ind = 0
        for sample in self.atomList.instance.positives:
            score = self.matching_score(sample)
            self.positive_histogram[score] += 1
            self.pos_error_indexes[score].append(ind)
            ind += 1

        # mathcing score on negative samples
        ind = 0
        for sample in self.atomList.instance.negatives:
            score = self.matching_score(sample)
            self.negative_histogram[score] += 1
            self.neg_error_indexes[score].append(ind)
            ind += 1
    #---------------------------------------------------------------------------
    def matching_score(self, sample):
        res = 0
        for ind_atom in self.body: # up to now, multiple atoms is not handled
            if sample[self.atomList.atoms[ind_atom][0]] != self.atomList.atoms[ind_atom][1]:
                res += 1
        return res

    #---------------------------------------------------------------------------
    def body_length(self):
        return len(self.body)

    #---------------------------------------------------------------------------
    def is_dominated(self, other):
        res = False
        if self.score[0] == other.score[0] and self.score[1] == other.score[1]:
            res = False
        else:
            if self.score[0] >= other.score[0] and self.score[1] <= other.score[1]:
                res = True
        return res
    #---------------------------------------------------------------------------
    def to_string(self):
        res = ''
        for ind in range(len(self.body)):
            res += str(self.atomList.atoms[self.body[ind]][0])
            res += '_'
            res += str(self.atomList.atoms[self.body[ind]][1])
            if ind < len(self.body)-1:
                res += ', '
        return res

    #---------------------------------------------------------------------------
    def histograms_to_string(self):
        res = ''

        barLength = 70

        res += "Negative histogram:\n"
        for ind_error in range(len(self.negative_histogram)):
            if ind_error < 10:
                res += '0'
            res += str(ind_error) + ': '
            freq = self.negative_histogram[ind_error]/len(self.atomList.instance.negatives)
            for step in range(math.floor(freq*barLength)):
                res += '#'
            res += '\n'

        res += "Positive histogram:\n"
        for ind_error in range(len(self.positive_histogram)):
            if ind_error < 10:
                res += '0'
            res += str(ind_error) + ': '
            freq = self.positive_histogram[ind_error]/len(self.atomList.instance.positives)
            for step in range(math.floor(freq*barLength)):
                res += '#'
            res += '\n'

        return res
