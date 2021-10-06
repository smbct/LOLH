#!/usr/bin/python

import numpy as np

#-------------------------------------------------------------------------------
def compute_subsets(elements, size):

    result = []
    stack = [0]
    backtrack = False

    while len(stack) > 0:

        # backtracking: go back to the last element that can still be used
        while backtrack:
            if len(stack) == 0:
                backtrack = False
            else:
                if stack[-1] < len(elements)-1 and len(elements)-stack[-1]-1 > size-len(stack):
                    backtrack = False
                    stack[-1] += 1
                else:
                    stack.pop()

        if len(stack) == size:
            # result.append(stack.copy())
            result.append([elements[ind] for ind in stack])
            backtrack = True

        elif len(stack) > 0:
            ind = stack[-1]
            if len(elements)-ind >= size-len(stack):
                stack.append(ind+1)
        else:
            backtrack = True

    return result

#-------------------------------------------------------------------------------
def compute_subsets_count(elements, size):

    # compute all subsets of a given size from a list of occurences
    # -> assume all values are unique and represented by their index
    # elements: (nbelt0, nbelt1, ...)

    result = []
    backtrack = False

    # number of itmes that are still available
    n_remaining = [0 for _ in elements]
    for ind in range(len(elements)-2, -1, -1):
        n_remaining[ind] = n_remaining[ind+1]+elements[ind+1]

    # start with index = 0
    ind = 0
    new_count = max(1,size-n_remaining[ind])
    stack = [(ind,new_count)]
    n_selected = new_count # number of items that have been selectd

    # cpt = 0

    while len(stack) > 0:

        # cpt += 1

        # print('current stack: ' + str(stack))
        # print('selected: ' + str(n_selected))
        # print('backtrack ? ' + str(backtrack))

        if backtrack:
            while backtrack:
                if len(stack) > 0:
                    ind,count = stack[-1]

                    if n_selected < size and count < elements[ind]: # increase the counter of the last element
                        stack[-1] = (ind,count+1)
                        n_selected += 1
                        backtrack = False
                    else: # otherwise pop the last element

                        stack.pop()
                        n_selected -= count

                        if ind < len(elements) - 1 and n_remaining[ind] >= size-n_selected: # add the next one if possible
                            new_count = max(1, size-n_selected-n_remaining[ind+1])
                            n_selected += new_count
                            stack.append((ind+1,new_count)) # take into account the remaining number of items for next initialization
                            backtrack = False
                else:
                    backtrack = False

        elif n_selected == size:
            result.append(stack.copy())
            backtrack = True
        else: # n_selected < size -> exploration

            ind,count = stack[-1]

            if ind < len(elements)-1 and n_remaining[ind] >= size-n_selected:
                new_count = max(1, size-n_selected-n_remaining[ind+1])
                stack.append((ind+1, new_count))
                n_selected += new_count
            else:
                backtrack = True

    # print('nb iterations: ' + str(cpt))

    return result


#-------------------------------------------------------------------------------
def create_association_list(elements, values):
    # create a list that associates each value in values to a set of elements with this value

    unique_values = []
    associated_elt = []

    for input_index in range(len(values)):
        # look for values[index] in result
        output_index = 0
        while output_index < len(unique_values) and not np.array_equal(unique_values[output_index],values[input_index]):
            output_index += 1

        if output_index < len(unique_values):
            associated_elt[output_index].append(elements[input_index])
        else: # add a new entry to the output
            unique_values.append(values[input_index])
            associated_elt.append([elements[input_index]])

    return unique_values, associated_elt


#-------------------------------------------------------------------------------
def compute_const_int(indexes, values, position):

    # compute an interval of constant values in elements such that index is in the interval
    # indexes: list of indexes from values list
    left_ind = position
    right_ind = position

    while (left_ind-1 >= 0) and (values[indexes[left_ind-1]] == values[indexes[position]]):
        left_ind -= 1

    while right_ind+1 < len(indexes) and values[indexes[right_ind+1]] == values[indexes[position]]:
        right_ind += 1

    return (left_ind, right_ind)



###################################################
# test compute
# my_list = [1,2,3,4,5,6,7,8,9,10]
# my_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#
# # print(my_list)
#
# result = compute_subsets(my_list, 3)
# print(result)
#
# for elt in result:
#     bin_vec = [0 for _ in range(len(my_list))]
#     for ind in elt:
#         bin_vec[ind] = 1
#     print(bin_vec)



###################################################
# test compution constant interval
#
# indexes = [0,1,2,3,4,5,6,7,8,9]
# values = [1,2,3,4,9,9,5,6,7]
#
# print(compute_const_int(indexes, values, 3))



##################################################
# test compute subset count

# my_list = [5, 10, 4, 6]
#
# print('element counts:')
# print(my_list)
#
# res = compute_subsets_count(my_list, 23)
#
# print('\nresult:')
# for elt in res:
#     print(elt)


##################################################
# test association list
#
# elements = [1,2,3,4,5,6,7]
# values =   [3,3,1,4,5,5,3]
#
# print(create_association_list(elements, values))
