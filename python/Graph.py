#!/usr/bin/python

import numpy as np

# Manage a neibourhood graph from a single cell dataset
#-------------------------------------------------------------------------------
class Graph:

    #---------------------------------------------------------------------------
    def __init__(self, file_name=None):

        # list of barcodes string
        self.barcodes = []

        # dictionary of indexes given barcode string
        self.barcode_inds = {}

        # list of successors and weights
        self.suc_inds = []
        self.suc_weights = []

        if not file_name is None:
            self.load_from_file(file_name)

        return

    #---------------------------------------------------------------------------
    def load_from_file(self, file_name):
        file = open(file_name, 'r')
        content = file.read()
        file.close()
        content = content.splitlines()
        for line in content:
            line = line.split(',')
            pred_ind = self.get_ind(line[0])
            line = line[1:]
            for elt in line:
                elt = elt.split(':')
                suc_ind = self.get_ind(elt[0])
                self.suc_inds[pred_ind].append(suc_ind)
                self.suc_weights[pred_ind].append(float(elt[1]))
        return

    #---------------------------------------------------------------------------
    def write_to_file(self, file_name):
        file = open(file_name, 'w')
        for pred_ind in range(len(self.suc_inds)):
            file.write(self.barcodes[pred_ind])
            for suc_ind in self.suc_inds[pred_ind]:
                file.write(',')
                file.write(self.barcodes[suc_ind])
                file.write(':')
                file.write(str(self.suc_weights[pred_ind]))
            file.write('\n')
        file.close()

    #---------------------------------------------------------------------------
    def get_ind(self, barcode):
        if barcode not in self.barcode_inds:
            self.barcode_inds[barcode] = len(self.barcodes)
            self.barcodes.append(barcode)
            self.suc_inds.append([])
            self.suc_weights.append([])
        return self.barcode_inds[barcode]


    #---------------------------------------------------------------------------
    def get_successors(self, barcode):
        return [self.barcodes[suc_ind] for suc_ind in self.suc_inds[self.barcode_inds[barcode]]]

    #---------------------------------------------------------------------------
    def get_successor_weights(self, barcode):
        return [suc_weight for suc_weight in self.suc_weights[self.barcode_inds[barcode]]]

    #---------------------------------------------------------------------------
    def build_from_dataframe(self, dataframe):
        return

    #---------------------------------------------------------------------------
    def to_string(self):
        res = ''
        res += 'n vertices: ' + str(len(self.barcodes))
        return res
