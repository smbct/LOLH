#!/usr/bin/python3

import itertools
import numpy as np
import math

from instance import Instance
import histogram
import utils

# linear programming solver: compute only supported solutions of bi-objective un-constrained problems
# Main algorithm based on a sort


class Solver:

    """
    A class implementing functions to optimize classification rules

    ...

    Attributes
    ----------
    instance : Instance
        instance of the classification problem

    Methods
    -------

    """

    #---------------------------------------------------------------------------
    def __init__(self, instance):

        self.instance = instance

        return

    #---------------------------------------------------------------------------
    def solve_extreme(self, body_length, objective, return_value = 0):

        # compute the optimal rules according to one of the score

        # return the score of the optimal solution on that objective: solution not dominated

        # if required, return the body, can be under three forms:
        #   - 0: no body returned
        #   - 1: one list of atoms that have this optimal score
        #   - 2: one list of atoms that are in the body plus a list of atoms with the same value that complete the body
        #   - 3: the list of all possible solution that have the same optimal value

        # create the comparison function
        comp_func = None

        if objective == 0: # optimize on the positive score
            comp_func = lambda x: (self.instance.atom_score[x][0], self.instance.n_negatives()-self.instance.atom_score[x][1])
        else: # optimize on the negative score
            comp_func = lambda x: (self.instance.n_negatives()-self.instance.atom_score[x][1], self.instance.atom_score[x][0])

        # sort the atoms to obtain the optimal solution
        sorted_atoms = self.select_unique_val(self.instance.atom_score, objective+1)
        # sorted_atoms = [index for index in range(self.instance.n_atoms())]
        sorted_atoms.sort(key = comp_func)

        # print([self.instance.atom_score[ind] for ind in range(2*body_length)])

        # print([self.instance.atom_score[ind] for ind in sorted_atoms[:2*body_length]])

        # create a solution based on the sorted list of atoms and on the score

        # compute the multi-objective value
        output_point = np.sum([self.instance.atom_score[index] for index in sorted_atoms[:body_length]], axis=0)

        output_bodies = None # c-style declaration

        if return_value == 1: # return one possible body
            output_bodies = [atom_index for atom_index in sorted_atoms[:body_length]]

        elif return_value > 1:

            # if there are tie at the end, this lead to several equivalent rules
            left_ind, right_ind = utils.compute_const_int(sorted_atoms, self.instance.atom_score, body_length-1)

            base_body = None
            choice_part = None

            if left_ind == right_ind:
                base_body = [atom_index for atom_index in sorted_atom[:body_length]] # part of the body identical for all rules
                choice_part = []
            else:
                base_body = [atom_index for atom_index in sorted_atoms[:left_ind]] # part of the body identical for all rules
                choice_part = [atom_index for atom_index in sorted_atoms[left_ind:right_ind+1]]

            if return_value == 2:
                output_bodies = (base_body, choice_part)
            elif return_value == 3:
                subsets = utils.compute_subsets(choice_part, body_length-len(base_body))
                output_bodies = [base_body + elt for elt in subsets]

        if return_value == 0:
            return output_point
        else:
            return (output_point, output_bodies)


    #---------------------------------------------------------------------------
    def solve_scalarization(self, body_length, lbda1, lbda2, return_body = 0):

        # return a list of pareto-optimal points for the given scalazitation
        # the rules associated to the value computed can be returned in different forms
        # complete list of rules or list of atoms that can be assembled to obtain the rules

        # return value:
        #    - 0: list of bi-objective values only
        #    - 1: list of values + possible bodies: list of indices with number of element to select in each list
        #    - 2: list of values + all possible bodies for these values
        #    - 3: one body for each point computed

        # May return duplicates in the points

        point_list = [] # list of points in the bi-objective space to return
        body_list = None # list of rule bodies returned, if specified
        if return_body > 0:
            body_list = []

        # solve a scalarization of the two objectives : max lambda1 * obj1 + labmda2 * obj2
        scores = [ lbda1 * self.instance.atom_score[ind][0] + lbda2 * self.instance.atom_score[ind][1] for ind in range(self.instance.n_atoms())]

        # keep only one atom per variable: make sure there is not two variables with same body
        # sorted_atoms = [index for index in range(self.instance.n_atoms())]
        sorted_atoms = self.select_unique_val(scores)

        sorted_atoms.sort(key = lambda x: scores[x], reverse=True)

        # print([ (ind,scores[ind]) for ind in sorted_atoms[:body_length]])
        # print([ (scores[ind]-scores[ind+1]) for ind in sorted_atoms[:100]]) # display the scores, beyond the body length

        # compute the tie to obtain all incomparable solutions
        left_ind, right_ind = utils.compute_const_int(sorted_atoms, scores, body_length-1)

        # print('equiv: ', (left_ind, right_ind))

        # compute the list of all objective values that are optimal for the scalarization

        if right_ind > body_length-1: # check if it is possible to obtain more than one solution

            # atoms part of the all the bodies
            base_body = [atom_index for atom_index in sorted_atoms[:left_ind]]

            partial_score = None # compute partial objective value from the comparable items
            if left_ind > 0:
                partial_score = np.sum([self.instance.atom_score[index] for index in base_body], axis=0)
            else:
                partial_score = [0,0]

            # list of atoms with the same scalarized value
            equivalent_atoms = [atom_index for atom_index in sorted_atoms[left_ind:right_ind+1]]

            # print('base body: ', base_body, ' -> ', [scores[elt] for elt in base_body])
            # print('equivalent atoms: ', equivalent_atoms, ' -> ', [scores[elt] for elt in equivalent_atoms])

            # create an association list between equivalent atoms and their objective value
            unique, elements = utils.create_association_list(equivalent_atoms, [self.instance.atom_score[elt] for elt in equivalent_atoms])

            # sort points: from the minimum to the maximum on positive score
            sorted_point_indexes = [index for index in range(len(unique))]
            sorted_point_indexes.sort(key = lambda index: unique[index][0]) # sort indexes using first objective (positive score)
            unique = [unique[index] for index in sorted_point_indexes] # order the lists
            elements = [elements[index] for index in sorted_point_indexes]

            # create a list of occurences for each unique element
            counts = [len(elt) for elt in elements]

            # compute all subset of given size of elements in unique
            # print('nb elem: ', body_length-left_ind)
            subsets = utils.compute_subsets_count(counts, body_length-left_ind)

            # print('unique: ', unique)
            # print('counts: ', counts)
            #
            # print('elements: ', elements)
            # print('subsets: ', subsets)
            # print([(unique[elt[0]], elt[1], elements[elt[0]]) for subset in subsets for elt in subset])
            #
            # print('partial score: ', partial_score)

            # compute all non dominated points
            for subset in subsets:

                # compute the value of the solution in the bi-objective space
                point = np.sum([ np.multiply(unique[elt[0]], elt[1]) for elt in subset], axis=0)
                point = np.sum([partial_score, point], axis=0)

                point_list.append(point) # new bi-objective point is added to the list

                # print('body returned: ', [(elt[1], elements[elt[0]]) for elt in subset])

                if return_body == 2: # compute all the bodies for this objective value
                    sub_body = []
                    for elt in subset:
                        sub_body.append(utils.compute_subsets(elements[elt[0]], elt[1]))
                    bodies = []
                    sub_bodies = itertools.product(*sub_body)
                    for elt in sub_bodies:
                        new_body = base_body.copy()
                        for sub_elt in elt:
                            new_body += sub_elt
                        bodies.append(new_body)
                    body_list.append(bodies)

                elif return_body > 2: # compute only one body
                    body_selected = base_body.copy()
                    for elt in subset:
                        body_selected += elements[elt[0]][:elt[1]]
                    body_list.append(body_selected)

                elif return_body == 1: # add the bodies in two parts: base atoms and choice
                    body_list.append([(len(base_body), base_body)] + [(elt[1], elements[elt[0]]) for elt in subset])

        else: # there are no tie, only one body is computed

            # only one body and only one associated point
            point_list.append(np.sum([ self.instance.atom_score[atom_index] for atom_index in sorted_atoms[:body_length] ], axis=0))

            if return_body == 1: # return bi-objective value + body, use same structure as before
                body_list.append((body_length, [atom_index for atom_index in sorted_atoms[:body_length]]))
            elif return_body == 2:
                body_list.append([[atom_index for atom_index in sorted_atoms[:body_length]]])
            elif return_body == 3:
                body_list.append([atom_index for atom_index in sorted_atoms[:body_length]])

        if return_body == 0:
            return point_list
        else:
            return (point_list,body_list)


    #---------------------------------------------------------------------------
    def compute_supported(self, body_length, return_weight = 0, return_body = 0):

        # compute all supported solutions of the instance, for the given body length

        # several return values:
        # - point list: objective values
        # - weight list: list of weights to obtain each point
        # - body list: list of all bodies, in several form

        # compute first the optimal solutions on separate objectives
        optim_pos,body_pos = self.solve_extreme(body_length, 0, 1)
        optim_neg,body_neg = self.solve_extreme(body_length, 1, 1)

        # print('pos: ', optim_pos)
        # print('neg: ', optim_neg)
        #
        # print('pos scores: ', [self.instance.atom_score[ind] for ind in body_pos])

        weight_pos = (-1,0)
        weight_neg = (0,1)

        # print('optimum on individual objectives')
        # print(optim_pos)
        # print(optim_neg)

        # check if the problem is degenerated
        if np.array_equal(optim_pos, optim_neg):
            # print('degenerated problem')

            if return_body > 0 and return_weight > 0:
                return ([optim_pos], [body_pos], [weight_pos])
            elif return_body > 0:
                return ([optim_pos], [body_pos])
            elif return_weight > 0:
                return ([optim_pos], [weight_pos])
            else:
                return [optim_pos]

        else:

            # Note: weights here are just for homogeneization purpose, the algorithm assumes there are no weakly dominated points

            point_list = [optim_pos, optim_neg] # list of points in the bi-objective space
            weight_list = None
            body_list = None

            if return_weight > 0: # need to be sorted afterward
                weight_list = [(-1,0), (0,1)]

            if return_body > 0:
                body_list = [body_pos, body_neg]

            new_points = None
            new_bodies = None
            new_weights = None

            # a stack is used for the recursive call of the scalarization
            stack = [(optim_pos, optim_neg)]

            while len(stack) > 0:
                left, right = stack.pop()

                # compute the ponderations
                # new_obj = max lbda1 * pos_score + lbda2 * neg_score
                lbda1 = left[1]-right[1]
                lbda2 = right[0]-left[0]

                if return_body > 0:
                    new_points, new_bodies = self.solve_scalarization(body_length, lbda1, lbda2, 3)
                else:
                    new_points = self.solve_scalarization(body_length, lbda1, lbda2, 0)

                # new_points_temp is useless here as only the indexes are used
                new_points_temp,indexes = np.unique(new_points, return_index=True, axis=0)

                # filter out points that are already known
                filtered_indexes = [index for index in indexes if new_points[index][0] > left[0] and new_points[index][1] < right[1]]

                if len(filtered_indexes) > 0: # make sure new points have been found

                    new_points = [new_points[index] for index in filtered_indexes]

                    point_list = np.concatenate((point_list, new_points))

                    if return_body > 0:
                        new_bodies = [new_bodies[index] for index in filtered_indexes]
                        body_list = np.concatenate((body_list, new_bodies), axis=0)

                    if return_weight > 0:
                        new_weights = [(lbda1, lbda2) for _ in filtered_indexes]
                        weight_list = np.concatenate((weight_list, new_weights), axis=0)


                    pending = np.concatenate(([left], new_points, [right]))
                    for index in range(len(pending)-1):
                        stack.append((pending[index], pending[index+1]))

            # sort all the points according to the first component
            indexes = list(range(len(point_list)))
            indexes.sort(key=lambda index: point_list[index][0])

            point_list = [point_list[index] for index in indexes]

            if return_body > 0:
                body_list = [body_list[index] for index in indexes]

            if return_weight > 0: # assume it contains bi-objective values followed by the weights
                weight_list = [weight_list[index] for index in indexes]

            if return_weight > 0 and return_body > 0:
                return point_list, body_list, weight_list
            elif return_weight > 0:
                return point_list, weight_list
            elif return_body > 0:
                return point_list
                # should be return point_list, body_list
            else:
                return point_list


    #---------------------------------------------------------------------------
    def compute_target_cover(self, body_length, threshold, rate, return_weight = 0, return_body = 0):

        # compute the first supported that cover some proportion of positive example

        # compute first the optimal solutions on separate objectives
        optim_pos,body_pos = self.solve_extreme(body_length, 0, 1)
        optim_neg,body_neg = self.solve_extreme(body_length, 1, 1)

        covered_pos = histogram.Histogram(self.instance, body_pos, histogram.Histogram_type.POSITIVE).positive_covered(threshold)
        covered_neg = histogram.Histogram(self.instance, body_neg, histogram.Histogram_type.POSITIVE).positive_covered(threshold)

        if np.array_equal(optim_pos, optim_neg) or covered_pos < rate:

            if return_body > 0 and return_weight > 0:
                return (optim_pos, body_pos, weight_pos)
            elif return_body > 0:
                return (optim_pos, body_pos)
            elif return_weight > 0:
                return (optim_pos, weight_pos)
            else:
                return optim_pos

        else:

            point_list = [optim_pos, optim_neg] # list of points in the bi-objective space
            # number of positive samples covered for each body computed:
            covering_list = [covered_pos, covered_neg]

            if covering_list[0] < rate:
                print('cover rate: ', covering_list[0])
                print('cover rate neg: ', covering_list[1])

                raise('The covering rate cannot be acheived')

            weight_list = None
            body_list = None

            if return_weight > 0: # need to be sorted afterward
                weight_list = [(-1,0), (0,1)]

            if return_body > 0:
                body_list = [body_pos, body_neg]

            # keep track of all histograms computed so far
            positive_histograms = []

            # a stack is used for the recursive call of the scalarization
            # contains indexes of points
            stack = [(optim_pos, optim_neg)]

            # point indices instead of
            stack = [(0,1)]

            stop = False # true when the target point is found

            best_index = -1

            while not stop and len(stack) > 0:

                ind_left, ind_right = stack.pop()
                left, right = point_list[ind_left], point_list[ind_right]

                # compute the ponderations
                # new_obj = max lbda1 * pos_score + lbda2 * neg_score
                lbda1 = left[1]-right[1]
                lbda2 = right[0]-left[0]

                new_points, new_bodies = self.solve_scalarization(body_length, lbda1, lbda2, 3)

                # filter out points that are already known
                indexes = list(range(len(new_points)))
                filtered_indexes = [index for index in list(range(len(new_points))) if new_points[index][0] > left[0] and new_points[index][1] < right[1]]

                # sort the indexes
                filtered_indexes.sort(key = lambda index: new_points[index][0])

                new_points = [new_points[index] for index in filtered_indexes]
                new_bodies = [new_bodies[index] for index in filtered_indexes]
                new_histograms = [histogram.Histogram(self.instance, body, histogram.Histogram_type.POSITIVE) for body in new_bodies]
                new_weigths = [(lbda1, lbda2) for _ in range(len(new_points))]

                # look for the max and min points such that the target is inside (if there is a target)

                # compute positive samples covered by the body


                if len(new_points) > 0:

                    # create a list of points to compute the new weights
                    pending = np.concatenate(([ind_left], [len(point_list) + ind for ind in range(len(new_points))], [ind_right]))

                    # store the points with the weights -> to easily compute the body later if necessary

                    # add the points computed to the list of points already processed
                    point_list = np.concatenate((point_list, new_points))

                    if return_body > 0:
                        body_list = np.concatenate((body_list, new_bodies), axis=0)

                    covering_list = np.concatenate((covering_list, [histo.positive_covered(threshold) for histo in new_histograms]))
                    del new_histograms

                    if return_weight > 0:
                        weight_list = np.concatenate((weight_list, new_weigths), axis=0)

                    #if return_body > 0:
                    #    new_points = np.unique(new_points, axis=0)

                    if len(new_points) > 0: # recursive step: solve new scalarizations

                        # pending = np.concatenate(([left], new_points, [right]))

                        pending_index = len(pending)-1 # record the index of the first point verifying the target
                        found = False

                        while not found and pending_index >= 0:
                            # check if the target is in the interval
                            if covering_list[pending[pending_index]] >= rate:
                                best_index = pending[pending_index]
                                found = True
                                # print('best cover: ', covering_list[pending[pending_index]])
                            else:
                                pending_index -= 1

                        if found and pending_index < len(pending)-1: # add a new pair of weights
                            stack.append((pending[pending_index], pending[pending_index+1]))
                        # else:
                            # print('not found: ', pending_index, ' vs ', len(pending))
                            # print('left covering; ')

                else:
                    # no solution for the scalarization
                    # in each iteration, ind_left already satisfies the covering
                    best_index = ind_left


            # sort all the points according to the first component
            # indexes = list(range(len(point_list)))
            # indexes.sort(key=lambda index: point_list[index][0])
            # point_list = [point_list[index] for index in indexes]

            if return_weight > 0:
                return point_list[best_index], body_list[best_index], weight_list[best_index]
            else:
                return point_list[best_index], body_list[best_index]


    #---------------------------------------------------------------------------
    def select_unique_val(self, scores, bi_obj=0):

        # select only one atom per variable
        greater = None

        if bi_obj == 1: # positive error
            greater = lambda elt1, elt2: scores[elt1][0] < scores[elt2][0] or (scores[elt1][0] == scores[elt2][0] and scores[elt1][1] > scores[elt2][1])
        elif bi_obj == 2: # negative error
            greater = lambda elt1, elt2: scores[elt1][1] > scores[elt2][1] or (scores[elt1][1] == scores[elt2][1] and scores[elt1][0] < scores[elt2][0])
        else:
            greater = lambda elt1, elt2: scores[elt1] > scores[elt2]

        is_selected = [False for _ in range(self.instance.n_atoms())]

        for feature in self.instance.prediction_features:
            max = -1
            for value in range(self.instance.n_values[feature]):
                atom_index = self.instance.atom_indexes[(feature, value)]

                if max < 0 or greater(atom_index, max):
                    max = atom_index

            is_selected[max] = True

        return_value = [index for index in range(self.instance.n_atoms()) if is_selected[index]]

        del is_selected

        return return_value

    #---------------------------------------------------------------------------
    def select_best_atoms(self):

        # fast method to compute the "best atoms" to classify the instance
        # goal: compute a co-expression network


        # create line equations from these points
        x1, y1 = (0., 0.)
        x2, y2 = (float(self.instance.n_positives()), float(self.instance.n_negatives()))
        a = 1.
        b = (x1-x2)/(y2-y1)
        c = -a*x1 - b*y1

        max_dist = b*y2+c

        # compute the list of distances to the diagonal
        atom_distance = [-1 for ind_atom in range(self.instance.n_atoms())]
        for ind_atom in range(self.instance.n_atoms()):
            score = self.instance.atom_score[ind_atom]
            atom_distance[ind_atom] = (a*score[0] + b*score[1] + c)/max_dist

        sorted_atoms = [ind for ind in range(self.instance.n_atoms())]
        sorted_atoms.sort(key=lambda ind:atom_distance[ind], reverse=True)

        #sorted_atoms = sorted_atoms[:n_atoms_selected]

        scores = [atom_distance[index] for index in sorted_atoms]

        return sorted_atoms, scores

    #---------------------------------------------------------------------------
    def select_best_atoms_fast(self):

        # faster method to compute the "best atoms" to classify the instance
        # goal: compute a co-expression network


        # compute the list of distances to the diagonal
        atom_distance = [-self.instance.n_negatives()*self.instance.atom_score[ind_atom][0] + self.instance.n_positives()*self.instance.atom_score[ind_atom][1] for ind_atom in range(self.instance.n_atoms())]

        max_score = self.instance.n_negatives()
        min_score = -self.instance.n_negatives()

        atom_distance = [ (value-min_score)/(max_score-min_score) for value in atom_distance]

        sorted_atoms = [ind for ind in range(self.instance.n_atoms())]
        sorted_atoms.sort(key=lambda ind:atom_distance[ind], reverse=True)

        scores = [atom_distance[index] for index in sorted_atoms]

        return sorted_atoms, scores

    #---------------------------------------------------------------------------
    def select_k_best_atoms(self, n_atoms_selected):
        sorted_atoms, scores = self.select_best_atoms()
        return sorted_atoms[:n_atoms_selected], scores[n_atoms_selected]

    #---------------------------------------------------------------------------
    def select_best_atoms_threshold(self, threshold):

        sorted_atoms, scores = self.select_best_atoms()

        for ind in range(len(scores)):
            if scores[ind] >= threshold:
                ind += 1
            else:
                break

        return sorted_atoms[:ind], scores[:ind]

    #---------------------------------------------------------------------------
    def compute_nondominated_atoms(self, atoms):

        # return a list of atoms that are not dominated (based on their individual score)
        
        inst = self.instance
        
        sorted_atoms = atoms.copy()

        sorted_atoms.sort(key = lambda atom_ind: (inst.atom_score[atom_ind][0], -inst.atom_score[atom_ind][1]))

        nondominated_atoms = [sorted_atoms[0]]
        for ind_atom in sorted_atoms[1:]:
            
            new_score = inst.atom_score[ind_atom]
            last_score = inst.atom_score[nondominated_atoms[-1]]
            
            dominated = True
            
            if new_score[0] == last_score[0] and new_score[1] == last_score[1]:
                dominated = False
            
            if new_score[0] > last_score[0] and new_score[1] > last_score[1]:
                dominated = False
                
            if not dominated:
                nondominated_atoms.append(ind_atom)
                
        return nondominated_atoms
        
    #---------------------------------------------------------------------------
    def compute_relative_atom_area(self, nondominated_atoms):

        # compute the relative area from a list of nondominated atoms        

        inst = self.instance       

        prec_score = inst.atom_score[nondominated_atoms[0]]
    
        # first square, from the bottom to the first point
        relative_area = prec_score[1]*(inst.n_positives()-prec_score[0])
        
        for ind_atom in nondominated_atoms[1:]:
            current_score = inst.atom_score[ind_atom]
            relative_area += (current_score[1]-prec_score[1])*(inst.n_positives()-current_score[0])
            prec_score = current_score

        relative_area /= (inst.n_positives()*inst.n_negatives())

        return relative_area

    #---------------------------------------------------------------------------
    @staticmethod
    def relative_supported_area(solutions) -> float:

        points = solutions[0]

        complete_area = points[-1][0] - points[0][0]
        complete_area *= points[-1][1] - points[0][1]

        dominated_area = 0

        for ind_point in range(1,len(points)):
            width = points[-1][0] - points[ind_point][0]
            height = points[ind_point][1] - points[ind_point-1][1]
            dominated_area += width*height

        return float(dominated_area)/float(complete_area)


    #---------------------------------------------------------------------------
    @staticmethod
    def compute_diagonal_distances(solutions) -> list:

        # return the distance of each point of the pareto front to the diagonal

        # create the diagonal equation a*x+b*y+c = 0
        p1 = solutions[0]
        p2 = solutions[-1]

        if abs(p2[1] - p1[1]) > abs(p2[0] - p1[0]):
            ratio = (p2[1] - p1[1])/(p2[0] - p1[0])
            b = 1.
            a = -ratio
        else:
            ratio = (p2[0] - p1[0])/(p2[1] - p1[1])
            a = 1.
            b = -ratio
        c = -a*p1[0] -b*p1[1]

        distances = [0 for _ in range(len(solutions))]

        # compute all the distances
        k = math.sqrt(a**2+b**2)
        max_dist = 0.
        projected_point = None
        for ind_point in range(1,len(solutions)-1):

            point = solutions[ind_point]

            p3 = [0,0]
            p3[0] = (b*(b*point[0] - a*point[1]) - a*c)/(a*a+b*b)
            p3[1] = (a*(-b*point[0] + a*point[1]) -b*c)/(a*a+b*b)

            # d1 = math.sqrt( (point[0] - p3[0])**2 + (point[1] - p3[1])**2 )
            distances[ind_point] = abs(a*point[0]+b*point[1]+c)/k

        return distances


    #---------------------------------------------------------------------------
    @staticmethod
    def max_min_dist(solutions) -> (float,float):

        # return the max of the min distance of a point of the pareto front to the diagonal

        # create the diagonal equation a*x+b*y+c = 0
        p1 = solutions[0]
        p2 = solutions[-1]

        if abs(p2[1] - p1[1]) > abs(p2[0] - p1[0]):
            ratio = (p2[1] - p1[1])/(p2[0] - p1[0])
            b = 1.
            a = -ratio
        else:
            ratio = (p2[0] - p1[0])/(p2[1] - p1[1])
            a = 1.
            b = -ratio
        c = -a*p1[0] -b*p1[1]

        # compute all the distances
        k = math.sqrt(a**2+b**2)
        max_dist = 0.
        projected_point = None
        for ind_point in range(1,len(solutions)-1):

            point = solutions[ind_point]

            p3 = [0,0]
            p3[0] = (b*(b*point[0] - a*point[1]) - a*c)/(a*a+b*b)
            p3[1] = (a*(-b*point[0] + a*point[1]) -b*c)/(a*a+b*b)

            # d1 = math.sqrt( (point[0] - p3[0])**2 + (point[1] - p3[1])**2 )
            dist = abs(a*point[0]+b*point[1]+c)/k

            if dist > max_dist:
                max_dist = dist
                projected_point = p3

        # utopia point to obtain a relative distance
        utopia = (p1[0], p2[1])
        distUtopia = abs(a*utopia[0]+b*utopia[1]+c)/k
        max_dist /= distUtopia

        # compute the relative position of the projected point on the diagonal
        # 0 <-> 1 : pos_samples <-> neg_samples
        relative_pos = math.sqrt( (projected_point[0] - solutions[0][0])**2 + (projected_point[1] - solutions[0][1])**2)
        relative_pos = relative_pos / math.sqrt( (solutions[-1][0] - solutions[0][0])**2 + (solutions[-1][1] - solutions[0][1])**2)

        return max_dist, relative_pos


    #---------------------------------------------------------------------------
    @staticmethod
    def max_min_dist_ind(solutions) -> int:

        # same but return an index

        # create the diagonal equation a*x+b*y+c = 0
        p1 = solutions[0]
        p2 = solutions[-1]

        if abs(p2[1] - p1[1]) > abs(p2[0] - p1[0]):
            ratio = (p2[1] - p1[1])/(p2[0] - p1[0])
            b = 1.
            a = -ratio
        else:
            ratio = (p2[0] - p1[0])/(p2[1] - p1[1])
            a = 1.
            b = -ratio
        c = -a*p1[0] -b*p1[1]

        # compute all the distances
        k = math.sqrt(a**2+b**2)
        max_dist = 0.
        max_ind = 0
        projected_point = None
        for ind_point in range(1,len(solutions)-1):

            point = solutions[ind_point]

            p3 = [0,0]
            p3[0] = (b*(b*point[0] - a*point[1]) - a*c)/(a*a+b*b)
            p3[1] = (a*(-b*point[0] + a*point[1]) -b*c)/(a*a+b*b)

            # d1 = math.sqrt( (point[0] - p3[0])**2 + (point[1] - p3[1])**2 )
            dist = abs(a*point[0]+b*point[1]+c)/k

            if dist > max_dist:
                max_dist = dist
                max_ind = ind_point
                projected_point = p3


        return max_ind

    #---------------------------------------------------------------------------
    @staticmethod
    def mean_median(instance, points) -> (float,float):

        # return the media, of the positive mean score and the negative mean score

        # points are already sorted

        ind = int(math.floor(len(points)/2))

        if len(points) % 2 == 0:
            med_pos = (points[ind-1][0]/instance.n_positives()+points[ind][0]/instance.n_positives())/2
            med_neg = (points[ind-1][1]/instance.n_negatives()+points[ind][1]/instance.n_negatives())/2
        else:
            med_pos = points[ind][0]/instance.n_positives()
            med_neg = points[ind][1]/instance.n_negatives()

        return med_pos, med_neg


    #---------------------------------------------------------------------------
    @staticmethod
    def mean_var(instance, bodies, points) -> (float,float,float,float):

        mean_pos_mean = 0
        mean_pos_var = 0
        mean_neg_mean = 0
        mean_neg_var = 0

        # compute the mean and
        for ind in range(len(bodies)):

            # create the histogram
            histo = histogram.Histogram(instance, bodies[ind])

            # compute the mean
            pos_mean = points[ind][0]/instance.n_positives()
            neg_mean = points[ind][1]/instance.n_negatives()

            # alternative comp mean
            # pos_mean_2 = 0
            # neg_mean_2 = 0
            # for error in range(len(histo.positive_histogram)):
            #     pos_mean_2 += error*len(histo.positive_histogram[error])
            #     neg_mean_2 += error*len(histo.negative_histogram[error])
            #
            # pos_mean_2 /= instance.n_positives()
            # neg_mean_2 /= instance.n_negatives()
            #
            # print('mean comp: ')
            # print(pos_mean, ' ; ', pos_mean_2)
            # print(neg_mean, ' ; ', neg_mean_2)

            # compute tha variance
            pos_var = 0
            neg_var = 0
            for error in range(len(histo.positive_histogram)):
                pos_var += len(histo.positive_histogram[error])*(error-pos_mean)**2
                neg_var += len(histo.negative_histogram[error])*(error-neg_mean)**2
            pos_var /= instance.n_positives()
            neg_var /= instance.n_negatives()

            mean_pos_mean += pos_mean
            mean_pos_var += pos_var

            mean_neg_mean += neg_mean
            mean_neg_var += neg_var

        mean_pos_mean /= len(bodies)
        mean_pos_var /= len(bodies)
        mean_neg_mean /= len(bodies)
        mean_neg_var /= len(bodies)

        return mean_pos_mean, mean_pos_var, mean_neg_mean, mean_neg_var


    #---------------------------------------------------------------------------
    @staticmethod
    def mean_var_med(instance, bodies, points) -> (float,float,float,float):

        # same as before but compute medians instead
        pos_means = []
        pos_vars = []
        neg_means = []
        neg_vars = []

        # compute the mean and
        for ind in range(len(bodies)):

            # create the histogram
            histo = histogram.Histogram(instance, bodies[ind])

            # compute the mean
            pos_mean = points[ind][0]/instance.n_positives()
            neg_mean = points[ind][1]/instance.n_negatives()

            # compute tha variance
            pos_var = 0
            neg_var = 0
            for error in range(len(histo.positive_histogram)):
                pos_var += len(histo.positive_histogram[error])*(error-pos_mean)**2
                neg_var += len(histo.negative_histogram[error])*(error-neg_mean)**2
            pos_var /= instance.n_positives()
            neg_var /= instance.n_negatives()

            pos_means.append(pos_mean)
            pos_vars.append(pos_var)

            neg_means.append(neg_mean)
            neg_vars.append(neg_var)

        pos_means.sort()
        pos_vars.sort()
        neg_means.sort()
        neg_vars.sort()

        ind = int(math.floor(len(bodies)/2))

        if len(bodies)%2 == 0:
            med_pos_mean = (pos_means[ind-1]+pos_means[ind])/2.
            med_pos_var = (pos_vars[ind-1]+pos_vars[ind])/2.
            med_neg_mean = (neg_means[ind-1]+neg_means[ind])/2.
            med_neg_var = (neg_vars[ind-1]+neg_vars[ind])/2.
        else:
            med_pos_mean = pos_means[ind]
            med_pos_var = pos_vars[ind]
            med_neg_mean = neg_means[ind]
            med_neg_var = neg_vars[ind]

        return med_pos_mean, med_pos_var, med_neg_mean, med_neg_var

    #---------------------------------------------------------------------------
    @staticmethod
    def lift_instance(instance, body_length, pos_ratio, neg_ratio):

        # modify the instance to improve the pareto front

        pos_index = {}
        ind = 0
        for elt in instance._pos_samples:
            pos_index[elt] = ind
            ind += 1

        neg_index = {}
        ind = 0
        for elt in instance._neg_samples:
            neg_index[elt] = ind
            ind += 1

        if instance.n_positives() > 0 and instance.n_negatives() > 0:

            solver = Solver(instance)

            # compute supported solutions
            supported_solutions = solver.compute_supported(body_length, 1, 1)

            # compute relative area
            area = Solver.relative_supported_area(supported_solutions)
            print('(lifting) relative area: ', area)

            # compute sampl score
            pos_sample_error = [0 for _ in range(instance.n_positives())]
            neg_sample_error = [0 for _ in range(instance.n_negatives())]

            for body in supported_solutions[1]:

                # create the histogram
                histo = histogram.Histogram(instance, body)

                for error in range(len(histo.positive_histogram)):
                    for sample in histo.positive_histogram[error]:
                        pos_sample_error[pos_index[sample]] += error
                    for sample in histo.negative_histogram[error]:
                        neg_sample_error[neg_index[sample]] += error

            # sort the samples according to the cumulated scores
            sorted_pos = sorted(instance._pos_samples, key = lambda x: pos_sample_error[pos_index[x]], reverse=True)
            sorted_neg = sorted(instance._neg_samples, key = lambda x: neg_sample_error[neg_index[x]])

            ind_cut_pos = int(instance.n_positives()*pos_ratio)
            ind_cut_neg = int(instance.n_negatives()*neg_ratio)

            # print([pos_sample_error[pos_index[elt]] for elt in sorted_pos[:int(len(sorted_pos)/6)]])
            # print([neg_sample_error[neg_index[elt]] for elt in sorted_neg[:int(len(sorted_neg)/6)]])

        return sorted_pos[ind_cut_pos+1:] + sorted_neg[:ind_cut_neg], sorted_neg[ind_cut_neg+1:]+sorted_pos[ind_cut_pos+1:]
