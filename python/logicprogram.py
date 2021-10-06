#!/usr/bin/python

import itertools

import pandas as pd

import random

# Logic program data structure
# Manage logic rules, matching, generating data, etc

#-------------------------------------------------------------------------------
class Logic_Rule:

    #---------------------------------------------------------------------------
    def __init__(self, head, body):

        self.head = head
        self.body = body

        return

    #---------------------------------------------------------------------------
    def match(self, state) -> bool:

        # returns trud iif the body of the rule mathces the input state

        res = True

        ind_body = 0
        while res and ind_body < len(self.body):

            atom = self.body[ind_body]
            if state[atom[0]] != atom[1]:
                res = False
            else:
                ind_body += 1

        return res

    #---------------------------------------------------------------------------
    def to_string(self) -> str:
        res = str(self.head[0]) + ':ind_' + str(self.head[1]) + ' :-'
        for atom in self.body:
            res += ' ' + str(atom[0]) + ':ind_' + str(atom[1])
        return res

#-------------------------------------------------------------------------------
class Logic_Program:

    #---------------------------------------------------------------------------
    def __init__(self):

        self._variables = [] # labels for each variable
        self.variable_index = {}
        self._values = [] # number of discrete values for each variables

        self._atoms = [] # list of atoms: variables and values
        self._atom_indexes = {} # association list atoms to indexes

        self.rules = []

        return

    #---------------------------------------------------------------------------
    def to_string(self) -> str:
        res = ''
        for rule in self.rules:
            res += self._variables[rule.head[0]] + '_' + str(rule.head[1])
            res += ' :-'
            for elt in rule.body:
                res += ' ' + self._variables[elt[0]] + '_' + str(elt[1])
            res += '.\n'
        return res

    #---------------------------------------------------------------------------
    def add_rule(self, rule):
        self.rules.append(rule)
        return

    #---------------------------------------------------------------------------
    def get_atom_index(self, atom) -> int:
        return self._atom_indexes[atom]

    #---------------------------------------------------------------------------
    def export_to_file(self, file_name):

        file = open(file_name, 'w')

        # write the variables
        for ind_var in range(len(self._variables)):
            line = 'VAR ' + self._variables[ind_var]
            for val in range(self._values[ind_var]):
                line += ' ' + str(val)
            line += '\n'
            file.write(line)

        file.write('\n')

        # write the rules
        for rule in self.rules:

            line = self._variables[rule.head[0]]
            line += '(' + str(rule.head[1]) + ',T) :- '

            for ind_atom in range(len(rule.body)):
                atom = rule.body[ind_atom]
                line += self._variables[atom[0]] + '(' + str(atom[1]) + ',T-1)'
                if ind_atom < len(rule.body)-1:
                    line += ', '

            line += '.\n'

            file.write(line)

        file.close()

        return

    #---------------------------------------------------------------------------
    def generate_all_transitions(self):

        # generate all possible transitions in synchronous semantics

        # generate the set of discrete values for all variables
        values = [ list(range(self._values[var_ind])) for var_ind in range(len(self._variables))]
        # print(values)

        # geenrate the set of all possible states (!exponential)
        states = list(itertools.product(*values))
        # print(states)

        transitions = []

        for state in states:

            next = [[] for _ in range(len(self._variables))]

            # look for all possible conclusions from this state
            for rule in self.rules:

                if rule.match(state):
                    next[rule.head[0]].append(rule.head[1])

            # make sure the state is not a dead end
            dead = False
            for val in next:
                if len(val) == 0:
                    dead = True

            # if dead:
                # print(':( no possible evolution)')

            if not dead: # if there are som successors, make the transitions
                successors = itertools.product(*next)
                for suc in successors:
                    transitions.append((state, suc))

        return transitions

    #---------------------------------------------------------------------------
    def generate_all_transitions_csv_synchrone(self, file_name):

        # generate all possible transitions in synchronous semantics
        # returns the transitions as csv

        # 2 files:  - states of the system observed
        #           - transitions between stats

        # generate the set of discrete values for all variables
        values = [ list(range(self._values[var_ind])) for var_ind in range(len(self._variables))]
        # print(values)

        # geenrate the set of all possible states (!exponential)
        states = list(itertools.product(*values))
        # print(states)

        state_indexes = {}
        for ind in range(len(states)):
            state_indexes[states[ind]] = ind

        transitions = []

        for ind_state in range(len(states)):

            state = states[ind_state]

            next = [[] for _ in range(len(self._variables))]

            # look for all possible conclusions from this state
            for rule in self.rules:

                if rule.match(state) and not rule.head[1] in next[rule.head[0]]:
                    next[rule.head[0]].append(rule.head[1])



            # make sure the state is not a dead end
            dead = False
            for val in next:
                if len(val) == 0:
                    dead = True

            # if dead:
                # print(':( no possible evolution)')

            if not dead: # if there are som successors, make the transitions
                successors = itertools.product(*next)
                for suc in successors:
                    transitions.append((ind_state, state_indexes[suc]))

        # print('transitions: ', transitions)

        states_labels = ['s' + str(ind) for ind in range(len(states))]

        # create the states dataframe
        states_df = pd.DataFrame(states, index = states_labels, columns = self._variables)

        # create the transitions dataframe
        # print(states_df)
        states_df.to_csv(file_name + '_states.csv')

        transitions_labels = ['tr_' + str(ind) for ind in range(len(transitions))]

        transitions_map = [[states_labels[elt[0]],states_labels[elt[1]]] for elt in transitions]

        transitions_df = pd.DataFrame(transitions_map, index = transitions_labels, columns = ['T-1', 'T'])
        transitions_df.to_csv(file_name + '_transitions.csv')

        # print(transitions_df)

        return

    #---------------------------------------------------------------------------
    def generate_all_transitions_csv_asynchrone(self, file_name):

        # generate all possible transitions in synchronous semantics
        # returns the transitions as csv

        # 2 files:  - states of the system observed
        #           - transitions between stats

        # generate the set of discrete values for all variables
        values = [ list(range(self._values[var_ind])) for var_ind in range(len(self._variables))]
        # print(values)

        # geenrate the set of all possible states (!exponential)
        states = list(itertools.product(*values))
        # print(states)

        state_indexes = {}
        for ind in range(len(states)):
            state_indexes[states[ind]] = ind

        # random.shuffle(states)
        # states = states[:int(len(states)/2)]


        transitions = []

        for ind_state in range(len(states)):

            state = states[ind_state]

            next = [[] for _ in range(len(self._variables))]

            # look for all possible conclusions from this state
            for rule in self.rules:

                if rule.match(state) and not rule.head[1] in next[rule.head[0]]:
                    next[rule.head[0]].append(rule.head[1])


            # make sure the state is not a dead end
            dead = False
            ind_var = 0
            for values in next:
                if len(values) == 0:
                    dead = True
                else: # generate successors
                    for value in values:
                        if value != states[ind_state][ind_var]:
                            suc = list(states[ind_state])
                            suc[ind_var] = value
                            suc = tuple(suc)
                            transitions.append((ind_state, state_indexes[suc]))
                ind_var += 1

        # print('transitions: ', transitions)

        # pick a subset of transitions
        random.shuffle(transitions)
        transitions = transitions[:2000]

        states_indexes = [elt[i] for i in [0,1] for elt in transitions]
        states_indexes = list(set(states_indexes))
        print(len(states_indexes))

        states_labels = ['s' + str(ind) for ind in range(len(states))]

        # select subset of states from picked transitions

        # create the states dataframe
        states_df = pd.DataFrame([states[ind] for ind in states_indexes], index = [states_labels[ind] for ind in states_indexes], columns = self._variables)

        # create the transitions dataframe
        # print(states_df)
        states_df.to_csv(file_name + '_states.csv')

        transitions_labels = ['tr_' + str(ind) for ind in range(len(transitions))]

        transitions_map = [[states_labels[elt[0]],states_labels[elt[1]]] for elt in transitions]

        transitions_df = pd.DataFrame(transitions_map, index = transitions_labels, columns = ['T-1', 'T'])
        transitions_df.to_csv(file_name + '_transitions.csv')

        # print(transitions_df)

        return

    #---------------------------------------------------------------------------
    @staticmethod
    def create_from_dataframe(dataframe):

        # create an empty program based on variables from a state dataframe

        lp = Logic_Program()

        lp._variables = dataframe.columns.copy()
        for index in range(len(lp._variables)):

            feature = lp._variables[index]

            lp.variable_index[feature] = index

            values = dataframe[feature].unique()
            # values.sort()

            lp._values.append(len(values))

            for val in values:
                lp._atoms.append((feature, val))
                lp._atom_indexes[(feature, val)] = len(lp._atoms)-1

        return lp

    #---------------------------------------------------------------------------
    @staticmethod
    def load_from_file(file_name):

        lp = Logic_Program()

        file = open(file_name, 'r')
        content = file.read()
        content = content.splitlines()
        file.close()

        # print(content)

        for line in content:

            tokens = line.split(' ')

            # read the variables of the program
            if tokens[0] == 'VAR':
                lp._variables.append(tokens[1])
                lp._values.append(len(tokens[2:]))
                lp.variable_index[lp._variables[-1]] = len(lp._variables)-1

                # red the values
                for value in range(lp._values[-1]):
                    lp._atoms.append((lp._variables[-1], value))
                    lp._atom_indexes[lp._atoms[-1]] = len(lp._atoms)-1

            elif tokens[0] != '': # read the rules

                head = Logic_Program.extract_atom(tokens[0])
                head = (lp.variable_index[head[0]], head[1])

                body = []

                # extract all body atoms
                for token in tokens[2:]:
                    atom = Logic_Program.extract_atom(token)
                    body.append((lp.variable_index[atom[0]], atom[1]))

                lp.rules.append(Logic_Rule(head, body))

        # print(lp.to_string())

        return lp

    #---------------------------------------------------------------------------
    @staticmethod
    def extract_atom(token):


        ind1 = token.index('(')
        ind2 = token.index(',')

        ind3 = token.index('T', ind1)
        ind4 = token.index(')')

        var_label = token[0:ind1]

        var_value = int(token[ind1+1:ind2])

        delay = 0
        if ind3+1 < ind4:
            delay = int(token[ind3+1:ind4])

        return (var_label, var_value, delay)
