#!/usr/bin/python

import numpy as np
import pandas as pd

class Instance:

    """
    A class used to represent an Instance of a classification problem

    ...

    Attributes
    ----------
    dataset : pandas dataframe
        a dataframe containing all the features for all the variables
    prediction_features : list of str
        list of all the features used for the prediction
    target_feature : str
        the target feature to predict
    atom_indexes : dictionary (int,int) -> int
        index of an atom given its value
    atoms : list of (int,int)
        list of atom values: variable index plus discrete value

    _positive_samples : list of int
        list of indexes of the positive samples in the datafram
    _negative_samples : list of int
        todo

    Methods
    -------
    compute_score()
        Compute the positive and negative score of each atom
    """

    #---------------------------------------------------------------------------
    def __init__(self):

        self.atom_indexes = None # dictionary of atom indexes
        self.atoms = None # list of atoms: variable index and discrete value
        self.atom_score = None # positive and negative score for each atom

        self.dataset = None
        self.prediction_features = None
        self.target_feature = None
        self.n_values = None # number of discrete values for each features

        self._pos_samples = None # positive samples of the instance
        self._neg_samples = None # negative samples of the instance

        return


    #---------------------------------------------------------------------------
    def _init_instance(self):

        # separate the instance into two subsets: positive and negative samples
        positive_dataset = self.dataset.loc[self.dataset[self.target_feature] == True]
        negative_dataset = self.dataset.loc[self.dataset[self.target_feature] == False]

        # list of sample labels for each class
        self._pos_samples = positive_dataset.index.copy()
        self._neg_samples = negative_dataset.index.copy()

        # size of each class
        n_positive = positive_dataset.shape[0]
        n_negative = negative_dataset.shape[0]

        self.n_values = {} # discrete values for each feature

        atom_index = 0

        # iterate over all prediction features of the dataset
        for feature in self.prediction_features:

            # count unique values
            pos_unique, pos_counts = np.unique(positive_dataset[feature], return_counts=True)
            pos_occ = dict(zip(pos_unique, pos_counts))

            neg_unique, neg_counts = np.unique(negative_dataset[feature], return_counts=True)
            neg_occ = dict(zip(neg_unique, neg_counts))

            # assemble and sort discrete values of the feature
            discrete_values = np.concatenate((pos_unique, neg_unique))
            discrete_values = np.sort(np.unique(discrete_values))

            self.n_values[feature] = len(discrete_values)

            # print('max val: ', np.max(discrete_values))

            for value in range(np.max(discrete_values)+1):

                self.atoms.append((feature, value))
                self.atom_indexes[self.atoms[-1]] = atom_index

                # compute the score:
                pos_score = n_positive
                if value in pos_unique:
                    pos_score -= pos_occ[value]

                # compute negative score
                neg_score = n_negative
                if value in neg_unique:
                    neg_score -= neg_occ[value]

                self.atom_score.append([pos_score, neg_score])

                atom_index += 1

            # create an atom and fill the score for each value
            # for value in discrete_values:
            #     self.atoms.append((feature, value))
            #     self.atom_indexes[self.atoms[-1]] = atom_index
            #
            #     # compute the score:
            #     pos_score = n_positive
            #     if value in pos_unique:
            #         pos_score -= pos_occ[value]
            #
            #     # compute negative score
            #     neg_score = n_negative
            #     if value in neg_unique:
            #         neg_score -= neg_occ[value]
            #
            #     self.atom_score.append([pos_score, neg_score])
            #
            #     atom_index += 1

        self.atom_score = np.array(self.atom_score) # convert to numpy array

        del positive_dataset
        del negative_dataset

        return


    #---------------------------------------------------------------------------
    def _init_sub_instance(self, pos_excluded, neg_excluded = []):

        # init the scores of a sub instance by modifying the original score
        # separate the istance into two subsets: positive and negative samples

        ex_pos_dataset = None
        ex_neg_dataset = None

        sub_pos = len(pos_excluded) > 0
        sub_neg = len(neg_excluded) > 0

        if sub_pos:
            ex_pos_dataset = self.dataset.loc[pos_excluded]
            self.atom_score[:,0] -= len(pos_excluded)

        if sub_neg:
            ex_neg_dataset = self.dataset.loc[neg_excluded]
            self.atom_score[:,1] -= len(neg_excluded)

        # iterate over all prediction features of the dataset
        for feature in self.prediction_features:

            if sub_pos: # count unique values excluded

                pos_unique, pos_counts = np.unique(ex_pos_dataset[feature], return_counts=True)
                pos_occ = dict(zip(pos_unique, pos_counts))

                for value in pos_unique: # add the error that should not be substracted
                    atom_ind = self.atom_indexes[(feature, value)]
                    self.atom_score[atom_ind][0] += pos_occ[value]

            if sub_neg:

                neg_unique, neg_counts = np.unique(ex_neg_dataset[feature], return_counts=True)
                neg_occ = dict(zip(neg_unique, neg_counts))

                for value in neg_unique: # add the error that should not be substracted
                    atom_ind = self.atom_indexes[(feature, value)]
                    self.atom_score[atom_ind][1] += neg_occ[value]

        if sub_pos:
            del ex_pos_dataset

        if sub_neg:
            del ex_neg_dataset


        return

    #---------------------------------------------------------------------------
    def n_positives(self) -> int:
        return len(self._pos_samples)

    #---------------------------------------------------------------------------
    def n_negatives(self) -> int:
        return len(self._neg_samples)

    #---------------------------------------------------------------------------
    def n_atoms(self) -> int:
        return len(self.atoms)

    #---------------------------------------------------------------------------
    def get_atom(self, atom_index):
        return self.atoms[atom_index]

    #---------------------------------------------------------------------------
    def has_atom(self, atom):
        return atom in self.atoms

    #---------------------------------------------------------------------------
    def get_atom_index(self, atom):
        return self.atom_indexes[atom]

    #---------------------------------------------------------------------------
    def to_string(self):
        res = ''

        return res

    #---------------------------------------------------------------------------
    def compute_rule_error(self, rule, index) -> float:

        sample = self.dataset.loc[index]

        error = 0
        for atom_index in rule:
            atom = self.atoms[atom_index]
            if sample[atom[0]] != atom[1]:
                error += 1

        return error


    #---------------------------------------------------------------------------
    @staticmethod
    def create_instance(dataset, prediction_features, target_feature):

        # create an instance from a dataframe

        res = Instance() # create the new Instance object

        # associate a name and index to each atom
        res.atom_indexes = {} # dictionary of atom indexes
        res.atoms = [] # list of atoms: variable index and discrete value
        res.atom_score = [] # positive and negative score for each atom

        res._pos_samples = [] # positive samples of the instance
        res._neg_samples = [] # negative samples of the instance

        # initialization of atom attributes
        res.dataset = dataset
        res.prediction_features = prediction_features
        res.target_feature = target_feature

        res._init_instance() # initiallization of atom scores

        return res


    #---------------------------------------------------------------------------
    @staticmethod
    def create_sub_instance(main_instance, pos_excluded, neg_excluded = []):

        res = Instance()

        # only references are used here
        res.atoms = main_instance.atoms
        res.atom_indexes = main_instance.atom_indexes

        res.dataset = main_instance.dataset
        res.prediction_features = main_instance.prediction_features
        res.target_feature = main_instance.target_feature

        # copy are made for scores and sample classes
        res.atom_score = main_instance.atom_score.copy()

        res._pos_samples = [elt for elt in main_instance._pos_samples if not elt in pos_excluded]
        res._neg_samples = [elt for elt in main_instance._neg_samples if not elt in neg_excluded]

        res.n_values = main_instance.n_values

        res._init_sub_instance(pos_excluded, neg_excluded)

        return res

    #---------------------------------------------------------------------------
    @staticmethod
    def create_instance_explicit(dataset, pos_samples, neg_samples):

        # create the instance by explicitly stating positive and negative samples

        dataset['Class'] = False

        # initialize the class feature: positive and negative cells
        dataset.loc[pos_samples,'Class'] = True

        # creation of the associated classification instance
        prediction_features = dataset.columns.drop(['Class'])
        instance = Instance.create_instance(dataset, prediction_features, 'Class')

        return instance

    #---------------------------------------------------------------------------
    @staticmethod
    def create_cluster_instance(dataframe, clusters, label):

        # Join the dataframe with the cell labels
        dataframe = dataframe.join(clusters)

        dataframe['Class'] = (dataframe['Label'] == label)

        features = dataframe.columns
        prediction_features = features.drop(['Class', 'Label'])
        target_feature = 'Class'

        instance = Instance.create_instance(dataframe, prediction_features, target_feature)

        return instance

    #---------------------------------------------------------------------------
    @staticmethod
    def create_coexpression_instance(dataset, gene_label, value):

        dataset['Class'] = (dataset[gene_label] == value)

        features = dataset.columns
        prediction_features = features.drop(['Class', gene_label]) # the gene column is removed as well

        instance = Instance.create_instance(dataset, prediction_features, 'Class')

        return instance

    #---------------------------------------------------------------------------
    @staticmethod
    def create_regulation_instance(dataset, transitions, gene_label, value, ratio, pred_neq = 0):

        # pred_neq == 0: any predecessors
        # pred_neq == 1: only predecessors with a different value than the target value
        # pred_neq == 2: only predecessors with the same value

        # print('gene label: ', gene_label, ' ; value: ', value)
        # print('ratio: ', ratio)
        # print('pred neq: ', pred_neq)
        # print('n transitions: ', transitions.shape[0])

        dataset['Class'] = False

        positive_barcodes = []

        for barcode in dataset.index:

            # filter out the predecessors based on the target value
            if pred_neq == 0 or (pred_neq == 1 and dataset[gene_label][barcode] != value) or (pred_neq == 2 and dataset[gene_label][barcode] == value):

                # compute successors
                successors = (transitions.loc[transitions['T-1'] == barcode])['T']

                # compute the number of successors verifying the target value
                pos_suc = dataset.loc[successors,gene_label] == value

                if pos_suc.shape[0] > 0:

                    if ratio <= -0.5:
                        positive_barcodes.append(barcode)
                    else:
                        occ = pos_suc.value_counts()
                        suc_ratio = float(occ[True]) / float(pos_suc.shape[0])
                        if suc_ratio >= ratio:
                            positive_barcodes.append(barcode)


        # initialize the class feature: positive and negative cells
        dataset.loc[positive_barcodes,'Class'] = True

        # creation of the associated classification instance
        prediction_features = dataset.columns.drop(['Class'])
        instance = Instance.create_instance(dataset, prediction_features, 'Class')

        return instance

    #---------------------------------------------------------------------------
    @staticmethod
    def create_regulation_instance_delay(dataset, transitions, gene_label, value, ratio, delay, pred_neq = 0):

        print('test delay')

        # compute a graph on the transitions
        row_index = {}
        ind = 0
        for barcode in dataset.index.values:
            row_index[barcode] = ind
            ind += 1

        row_label = dataset.index.values.copy()

        # compute a graph of successors
        graph = [[] for _ in range(len(row_label))]
        for ind in range(len(row_label)):
            pred = row_label[ind]
            suc_labels = transitions.loc[transitions['T-1'] == pred]['T']
            for elt in suc_labels:
                graph[ind].append(row_index[elt])

        successors_delay = [[] for _ in range(len(row_label))] # list of successors at given distance

        # for each cell, compute successors at given distance (delay)
        for ind in range(len(row_label)):
            stack = [(ind,0)]
            while len(stack) > 0:
                top = stack.pop()
                if top[1] == delay:
                    successors_delay[ind].append(top[0])
                elif top[1] < delay:
                    for suc_ind in graph[top[0]]:
                        stack.append((suc_ind, top[1]+1))

        # compute the class vector based on the ratio
        dataset['Class'] = False
        for ind in range(len(row_label)):

            if pred_neq == 0 or (pred_neq == 1 and dataset[ind][gene_label] != value) or (pred_neq == 2 and dataset[ind][gene_label] == value):

                n_suc = len(successors_delay[ind])

                # compute the successors with the correct value
                suc_dataframe = dataset.iloc[successors_delay[ind],:]
                suc_dataframe = suc_dataframe[suc_dataframe[gene_label] == value]
                # print(suc_dataframe)

                if n_suc > 0:
                    if ratio <= -0.5 or suc_dataframe.shape[0]/n_suc >= ratio:
                        dataset.loc[row_label[ind],'Class'] = True

        # creation of the associated classification instance
        prediction_features = dataset.columns.drop(['Class'])
        instance = Instance.create_instance(dataset, prediction_features, 'Class')

        return instance


    #---------------------------------------------------------------------------
    @staticmethod
    def create_random_instance(dataset, proportion):

        # create a random instance from the dataset, having a given proportion of positive samples

        # don't forget the copy, otherwise, the cell order in the dataset is modified
        # as a result, the instance select the k first cells as positive -> well conditioned instance
        indexes = dataset.index.values.copy()
        # indexes = dataset.index.values

        np.random.shuffle(indexes)

        n_pos = int(proportion*float(len(indexes)))

        dataset['Class'] = False

        dataset.loc[indexes[:n_pos], 'Class'] = True

        prediction_features = dataset.columns.drop(['Class'])
        instance = Instance.create_instance(dataset, prediction_features, 'Class')

        return instance
