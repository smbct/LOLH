#!/usr/bin/python3

import numpy as np
import pandas as pd

from instance import Instance
from solver import Solver
import histogram

from logicprogram import Logic_Program
from logicprogram import Logic_Rule

# functions that learn logic dynamical discrete programs
# interface between lp solver and logic program data structure and panda dataset


# inference functions to create a logic program from an atom
#-------------------------------------------------------------------------------
def infer_rules(instance, body_length, cover_threshold, nb_predictors):

    # compute number of samples covered by each rules (approximative)
    cover_rate = int(instance.n_positives()/nb_predictors)

    if cover_rate < 1:
        raise('Cannot make a model that learn no samples')

    bodies = [] # list of bodies created for this atom

    excluded_samples = [] # list of positive samples excluded for the re-optimization

    stop = False

    solver = Solver(instance)

    sub_instance = instance

    while not stop:

        # print('excluded: ', excluded_samples)

        solver = Solver(sub_instance)
        n_remaining = instance.n_positives() - len(excluded_samples)

        res = solver.compute_target_cover(body_length, cover_threshold, min(n_remaining,cover_rate), 0, 1)

        histo = histogram.Histogram(sub_instance, res[1], histogram.Histogram_type.POSITIVE)

        n_covered = histo.positive_covered(cover_threshold)

        # the rule is registered only of the covering is not null
        if n_covered > 0:
            bodies.append(res[1]) # new body is added to the list
            new_excluded = [elt for tab in histo.positive_histogram[:cover_threshold+1] for elt in tab]
            excluded_samples += new_excluded

            # make sure this is not the first iteration
            if n_remaining < instance.n_positives():
                del sub_instance

            if len(excluded_samples) == instance.n_positives():
                stop = True
            else:
                sub_instance = Instance.create_sub_instance(instance, excluded_samples)

        if n_covered < min(n_remaining,cover_rate): # goal is not reached, the body is too large
            if body_length > 0:
                body_length -= 1 # might be necessary to do some post processing
            else:
                stop = True

        del histo
        del solver

    return bodies

#-------------------------------------------------------------------------------
def create_dataset_class(states_df, transitions_df, feature, value) -> pd.DataFrame:

    # create the class column in the dataset corresponding to the dynamical feature learned

    # the transitions dataset contains transitions: pairs of labels from states dataframe
    # 2 col : T-1 and T

    instance_dataframe = states_df

    instance_dataframe['Class'] = False

    for tr_index in transitions_df.index:
        origin = transitions_df.loc[tr_index,'T-1']
        target = transitions_df.loc[tr_index,'T']
        if instance_dataframe.loc[target,feature] == value:
            instance_dataframe.loc[origin,'Class'] = True

    return instance_dataframe

#-------------------------------------------------------------------------------
def learn_program(states_df, transitions_df, body_length, cover_threshold, nb_predictors):

    # create an empty program from a state dataframe (no transition information yet)
    program = Logic_Program.create_from_dataframe(states_df)

    prediction_features = states_df.columns.copy()
    body_learned = []
    target_feature = 'Class'
    states_df['Class'] = False

    inst = None

    for col in prediction_features:

        values = states_df[col].unique()

        for value in values:

            rule_head = (program.variable_index[col], value)

            print('currently learning ' + col + '_' + str(value))

            instance_dataframe = create_dataset_class(states_df, transitions_df, col, value)
            inst = Instance.create_instance(instance_dataframe, prediction_features, target_feature)
            res = infer_rules(inst, body_length, cover_threshold, nb_predictors)
            body_learned.append(res)

            for body in res:
                rule_body = [ (program.variable_index[inst.atoms[ind][0]], inst.atoms[ind][1]) for ind in body ]
                rule = Logic_Rule(rule_head, rule_body)
                program.add_rule(rule)

            # print('inference result: ', res)
            print('inference result: ', [[inst.atoms[ind] for ind in body] for body in res ])
            del inst

    print('Program learned:')
    print(program.to_string())

    return program
