#!/usr/bin/python

import pandas as pd
import numpy as np
import math

from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.backend_bases import MouseButton # Interactive part

from instance import Instance
from solver import Solver
import histogram

import visualizer



#-------------------------------------------------------------------------------
# Global Variables
#-------------------------------------------------------------------------------

is_mouse_pressed = False
selection_coordinates = [None, None]
instance = None
supported_solutions = None
selection_rectangle = None
selection_indexes = []

instance_info = None

normalize = True
supported_points = []

#-------------------------------------------------------------------------------
# Callback functions

#-------------------------------------------------------------------------------
def onclick(event):

    global is_mouse_pressed
    global selection_coordinates

    if event.button == MouseButton.MIDDLE:
        is_mouse_pressed = True
        selection_coordinates[0] = [event.xdata, event.ydata]

    return

#-------------------------------------------------------------------------------
def onrelease(event):

    global is_mouse_pressed
    global selection_coordinates
    global selection_indexes

    global df_normalised
    global instance

    if event.button == MouseButton.MIDDLE:
        is_mouse_pressed = False
        selection_coordinates[1] = [event.xdata, event.ydata]
        update_indexes()
        redraw_pareto_plot()
        draw_histograms()

        # draw one violin plot (make sur the selection is not empty)
        if len(selection_indexes) > 0:
            fig,ax = plt.subplots()
            visualizer.plot_violins(ax, df_normalised, instance, supported_solutions[1][selection_indexes[0]])

            # scatter plots of individual gene scores
            fig,ax = plt.subplots()
            ax.set_xlabel('genes positive scores')
            ax.set_ylabel('genes negative scores')
            ax.set_title('Gene correlations w.r.t. NKG7_1')
            atom_colors = ['red' if ind in supported_solutions[1][selection_indexes[0]] else 'royalblue' for ind in range(instance.n_atoms())]
            ax.scatter([instance.atom_score[ind][0] for ind in range(instance.n_atoms())], [instance.atom_score[ind][1] for ind in range(instance.n_atoms())], marker='x', color=atom_colors)

            for ind in supported_solutions[1][selection_indexes[0]]:
                ax.text(instance.atom_score[ind][0] + 10, instance.atom_score[ind][1] - 10, instance.get_atom(ind)[0] + '_' + str(instance.get_atom(ind)[1]))


            plt.show(block=False)

    return

#-------------------------------------------------------------------------------
def onmove(event):

    if event.inaxes and is_mouse_pressed:
        ax = event.inaxes  # the axes instance
        selection_coordinates[1] = [event.xdata, event.ydata]
        redraw_pareto_plot()

    return

#-------------------------------------------------------------------------------
def redraw_pareto_plot():

    global is_mouse_pressed
    global supported_solutions
    global selection_coordinates
    global selection_rectangle
    global main_ax
    global main_fig

    global normalize
    global supported_points

    global instance_info

    if supported_solutions != None:

        main_ax.clear()

        if selection_coordinates[1] is None:
            visualizer.plot_rule_scores(main_ax, supported_points, None, instance_info)
        else:
            visualizer.plot_rule_scores(main_ax, supported_points, selection_indexes, instance_info)

        if not selection_coordinates[1] is None:

            coord = compute_rectangle()
            selection_rectangle = coord

            rect = patches.Rectangle((coord[0],coord[1]),coord[2],coord[3],linewidth=1,edgecolor='r',facecolor='none')
            main_ax.add_patch(rect)

        main_fig.canvas.draw()
        main_fig.canvas.flush_events()

    return

#-------------------------------------------------------------------------------
# compute bottom left coord and with/heigth of a rectangle
def compute_rectangle():
    global selection_coordinates
    left = min(selection_coordinates[0][0], selection_coordinates[1][0])
    bottom = min(selection_coordinates[0][1], selection_coordinates[1][1])
    width = abs(selection_coordinates[0][0]-selection_coordinates[1][0])
    height = abs(selection_coordinates[0][1]-selection_coordinates[1][1])
    return left,bottom,width,height

#-------------------------------------------------------------------------------
# Select the solutions corresponding
def update_indexes():

    global supported_solutions
    global selection_rectangle
    global selection_indexes

    global supported_points

    global instance

    selection_indexes = []

    for ind in range(len(supported_points)):
        point = supported_points[ind]
        if point[0] >= selection_rectangle[0] and point[1] >= selection_rectangle[1]:
            if point[0] <= selection_rectangle[0] + selection_rectangle[2]:
                if point[1] <= selection_rectangle[1] + selection_rectangle[3]:
                    selection_indexes.append(ind)

    # display the selected solutions
    print(selection_indexes)
    for ind in selection_indexes:
        print([instance.get_atom(ind2) for ind2 in supported_solutions[1][ind]])
    print('\n\n\n\n')

    return


#-------------------------------------------------------------------------------
# Plot the histograms corresponding to the selected solutions
def draw_histograms():

    global supported_solutions
    global selection_indexes
    global instance

    global instance_info

    n_histo = len(selection_indexes)

    width = int(math.sqrt(n_histo))
    if width*width < n_histo:
        width += 1

    if len(selection_indexes) > 0:

        fig, axs = plt.subplots(width, width, constrained_layout=True)

        axes = []
        if width > 1:
            axes = axs.flat[:n_histo]
        else:
            axes = [axs]

        ind_plot = 0
        for ax in axes:
            histo = histogram.Histogram(instance, supported_solutions[1][selection_indexes[ind_plot]])
            visualizer.plot_histograms(ax, histo, True, instance_info)
            ind_plot += 1

        plt.show(block=False)



#-------------------------------------------------------------------------------
# Connections between pyplot and callback functions
# Create the main plot
main_fig, main_ax = plt.subplots()

cid = main_fig.canvas.mpl_connect('button_press_event', onclick)
cid = main_fig.canvas.mpl_connect('button_release_event', onrelease)
cid = main_fig.canvas.mpl_connect('motion_notify_event', onmove)



#-------------------------------------------------------------------------------
# Main program

np.random.seed(42)



#-------------------------------------------------------------------------------
# MAIT cells
#-------------------------------------------------------------------------------


# # df = pd.read_csv('../dataset/MAIT/cells/binary_logcpm_reduced.csv', index_col=0) # Load the dataset: cell log count
# df = pd.read_csv('../dataset/MAIT/cells/binary_scanpy.csv', index_col=0) # scanpy cells
#
# df.rename_axis('Barcode', axis=0, inplace=True) # remove axis name (name of the previously first column)
# df_labels = pd.read_csv('../dataset/MAIT/cells_labels_corrected.csv', index_col=0, header=None, names=['Barcode', 'Label']) # Load the cluster label file
#
# # df_transitions = pd.read_csv('../dataset/MAIT/transitions_cells_ng.csv', index_col=0)
# # df_transitions = pd.read_csv('../dataset/MAIT/transitions/transitions_bidir_pt.csv', index_col=0) # bi-directional graph
#
# # df_transitions = pd.read_csv('../dataset/MAIT/transitions/transitions_n2_pt.csv', index_col=0) # pseudotime graph
# # df_transitions = pd.read_csv('../dataset/MAIT/transitions/transitions_n2_bidir.csv', index_col=0) # bi-directional graph
# # df_transitions = pd.read_csv('../dataset/MAIT/transitions/transitions_n2_rndpt.csv', index_col=0) # random pseudotime graph
# df_transitions = pd.read_csv('../dataset/MAIT/transitions/transitions_bidir_scanpy.csv', index_col=0) # scanpy
#
#
# instance = Instance.create_cluster_instance(df, df_labels, 'MAIT1')
# # instance = Instance.create_coexpression_instance(df, 'Gm45168', 1)
# # instance = Instance.create_regulation_instance(df, df_transitions, 'Shoc2', 0, 0.3, 1) # weird gene

#-------------------------------------------------------------------------------
# IMAGINE cells
#-------------------------------------------------------------------------------
file_name = '../dataset/IMAGINE/cell_types.csv'
df_cell_types = pd.read_csv(file_name, index_col=0)
df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)
print(df_cell_types.head())


file_name = '../../../Learning/IMAGINE_dataset/dataset/ctls_normalised_counts.csv'
df_normalised = pd.read_csv(file_name, index_col = 0).T
# print(df_normalised.head())

file_name = '../dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv'
# Load the dataset: cell log count
df = pd.read_csv(file_name, index_col=0)
# remove axis name (name of the previously first column)
df.rename_axis('Barcode', axis=0, inplace=True)
print(df.head())

#instance = Instance.create_cluster_instance(df, df_cell_types, 'CD8')

#-------------------------------------------------------------------------------
# instance = Instance.create_regulation_instance(df, df_transitions, 'Lmnb1', 1, 0.5, 1) # weird gene
# instance_info = ('Lmnb1', 1)

#-------------------------------------------------------------------------------
# instance = Instance.create_coexpression_instance(df, 'Gm45168', 1)
# instance = Instance.create_coexpression_instance(df, df.columns[0], 0) # index 0
# instance = Instance.create_coexpression_instance(df, 'Rorc', 1) # 1 gn d'1Traie
# instance = Instance.create_coexpression_instance(df, 'Ckb', 1) # 0.72

# instance = Instance.create_coexpression_instance(df, 'Lmo4', 1)

# instance = Instance.create_coexpression_instance(df, df.columns.values[366], 1)
# instance = Instance.create_coexpression_instance(df, 'Tyrobp', 1)

# print(df.columns.values[366])


#instance = Instance.create_coexpression_instance(df, 'Wdr76', 0) # 0.70
#instance_info = ('Rorc', 1)

instance = Instance.create_coexpression_instance(df, 'NKG7', 1)
instance_info = ('NKG7', 1)

#-------------------------------------------------------------------------------
# create random instance
# instance = Instance.create_random_instance(df, 0.1)
# instance_info = ('Rorc', 1)

print(instance_info)

body_length = 15

#--------------------------------------------------------
# Lifting
# print('Instance lifting')
# pos_samples, neg_samples = Solver.lift_instance(instance, body_length, 0.17, 0.17)
# instance_lifted = Instance.create_instance_explicit(df, pos_samples, neg_samples)
# instance = instance_lifted
#--------------------------------------------------------

print('Original instance:')
print('n atoms: ', instance.n_atoms())
print('n positive samples: ', instance.n_positives())
print('n negative samples: ', instance.n_negatives())

if instance.n_positives() > 0 and instance.n_negatives() > 0:

    solver = Solver(instance)

    supported_solutions = [ [], [], [] ]
    # for body_length in [5, 10, 20, 30, 40, 50, 60]:
    # for body_length in [10, 20, 30, 40]:
    for body_length in [15]:

        # compute supported solutions
        supported_solutions_temp = solver.compute_supported(body_length, 1, 1)
        supported_solutions[0] += supported_solutions_temp[0]
        supported_solutions[1] += supported_solutions_temp[1]
        supported_solutions[2] += supported_solutions_temp[2]

        selected_body = -1
        best_dist = -1
        for ind in range(len(supported_solutions_temp[0])):
            if selected_body == -1 or abs(supported_solutions_temp[0][ind][1] - supported_solutions_temp[0][ind][0]) > best_dist:
                best_dist = abs(supported_solutions_temp[0][ind][1] - supported_solutions_temp[0][ind][0])
                selected_body = ind

        print( [ instance.get_atom(ind) for ind in supported_solutions_temp[1][selected_body] ])

        # compute quality measures
        # area = Solver.relative_supported_area(supported_solutions_temp)
        # print('relative area: ', area)
        #
        # max_dist, dist_pos = Solver.max_min_dist(supported_solutions_temp[0])
        # print("max dist: ", max_dist)
        # print("max dist pos: ", dist_pos)
        #
        # pos_med, neg_med = Solver.mean_median(instance, supported_solutions_temp[0])
        # print('pos mean median: ', pos_med)
        # print('neg mean medan: ', neg_med)

    if normalize == True:
        supported_points = [ (supported_solutions[0][ind][0] / len(supported_solutions[1][ind]), supported_solutions[0][ind][1] / len(supported_solutions[1][ind])) for ind in range(len(supported_solutions[0]))]
    else:
        supported_points = supported_solutions[0]

    # compute relative area
    # area = Solver.relative_supported_area(supported_solutions)
    # print('relative area: ', area)

    main_ax.set_title('Supported points')
    redraw_pareto_plot()

    plt.show()

else:

    print('ERROR: null proportion of positive/negative samples')
