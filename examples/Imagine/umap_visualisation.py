#!/usr/bin/python3

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.patches import Patch


# LOLH objects
import sys
sys.path.append('../../python')

import solver

import visualizer
import histogram


#-------------------------------------------------------------------------------
def plot_cell_types():

    # import the cell types
    file_name = '../../dataset/Imagine/cell_types.csv'
    df_cell_types = pd.read_csv(file_name, index_col=0)
    df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)

    # import the UMAP 2d representation
    file_name = '../../dataset/Imagine/umap_coordinates.csv'
    umap_coordinates = pd.read_csv(file_name, index_col = 0)

    colors = [(248, 118, 109), (232, 133, 38), (211, 146, 0), (183, 159, 0), (147, 170, 0), (94, 179, 0)]
    colors += [(0, 186, 56), (0, 191, 116), (0, 193, 159), (0, 191, 196), (0, 185, 227), (0, 173, 250)]
    colors += [(97, 156, 255), (174, 135, 255), (219, 114, 251), (245, 100, 227), (255, 97, 195), (255, 105, 156)]

    for index in range(len(colors)):
        col = colors[index]
        colors[index] = (col[0]/255., col[1]/255., col[2]/255.)

    cell_types = ['CD4-naive', 'CD14', 'CD4', 'NK', 'CD8', 'CD16', 'B-naive', 'Treg']
    cell_types += ['gdT', 'CD8-naive', 'B', 'T_undef', 'cDC', 'platelet_cont', 'CD14-PPBP', 'HSC', 'pDC', 'undefined']

    cell_types_corrected = ['CD4-naive', 'CD14', 'CD4', 'NK', 'CD8', 'CD16', 'B-naive', 'T-regs']
    cell_types_corrected += ['T-gamma-delta', 'CD8-naive', 'B', 'undefined-T', 'cDC', 'contaminants', 'CD14-PPBP', 'HSC', 'pDC', 'undefined']

    # dictionary between cell type and color
    cell_color = {cell_types[ind]:colors[ind] for ind in range(len(colors))}

    # plot the UMAP with the cell types
    fig, ax = plt.subplots()
    ax.scatter(umap_coordinates['UMAP_1'][:], umap_coordinates['UMAP_2'][:], s=1, c = [cell_color[df_cell_types['Label'][cell]] for cell in umap_coordinates.index])
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    # ax.set_title('Visualization of the cell types')

    # print the cell types
    celltypes_index = 0
    for celltype in cell_types:
        if celltype != 'undefined':

            if celltype == 'T_undef':
                pos = (1.15, -0.15)
            elif celltype == 'Treg':
                pos = (-0.3, 1.3)
            elif celltype == 'CD4':
                pos = (-1.2, 3.6)
            elif celltype == 'CD8':
                pos = (5.5, -2.75)
            elif celltype == 'CD8-naive':
                pos = (-2, 9.15)
            elif celltype == 'B-naive':
                pos = (-14.7, 7.55)
            elif celltype == 'platelet_cont':
                pos = (-4.7, -2.2)
            else:
                pos = (0, 0)
                sub_cells = df_cell_types.loc[df_cell_types['Label'] == celltype].index
                for barcode in sub_cells:
                    pos = (pos[0] + umap_coordinates['UMAP_1'][barcode], pos[1]+umap_coordinates['UMAP_2'][barcode])
                pos = (pos[0] / len(sub_cells), pos[1] / len(sub_cells))
            ax.text(*pos, cell_types_corrected[celltypes_index], horizontalalignment='center', verticalalignment='center')

        celltypes_index += 1
    # fig.tight_layout(h_pad=1)

    legend_elements = [ Patch(facecolor=cell_color[cell_types[type_index]], label=cell_types_corrected[type_index]) for type_index in range(len(cell_types))]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., labelspacing=0.1)

    fig.tight_layout(h_pad=1)


    plt.show()

#-------------------------------------------------------------------------------
def plot_cell_macrotypes():

    # import the cell types
    file_name = '../../dataset/Imagine/cell_types_macro.csv'
    df_macro = pd.read_csv(file_name, index_col=0)
    df_macro.rename(columns={'cellType_macro': 'Label'}, inplace=True)

    macro_types = [str(label) for label in np.unique(list(df_macro['Label'].values))]
    macrotype_dic = {str(barcode):str(df_macro['Label'][barcode]) for barcode in df_macro.index}


    # import the UMAP 2d representation
    file_name = '../../dataset/Imagine/umap_coordinates.csv'
    umap_coordinates = pd.read_csv(file_name, index_col = 0)

     # define coordinates for label positions
    macro_bbox = {}
    macro_bbox['B'] = [(-16.5, -13.5), (2, 9)]
    macro_bbox['HSC'] = [(-6.5, -4.5), (6, 7)]
    macro_bbox['Myeloid'] = [(-8, 8), (-16, -11)]
    macro_bbox['T'] = [(-7.5, 14), (-7, 11)]
    macro_bbox['contamination'] = [(-6, -4), (-2.5, -1)]

    colors_dic = {}
    colors_dic['T'] = (0.973, 0.463, 0.427)
    colors_dic['Myeloid'] = (0.639, 0.647, 0.0)
    colors_dic['B'] = (0, 0.749, 0.49)
    colors_dic['contamination'] = (0, 0.69, 0.965)
    colors_dic['HSC'] = (0.906, 0.42, 0.953)
    colors_dic['nan'] = (0.3, 0.3, 0.3)

    fig, ax = plt.subplots()
    col = [colors_dic[macrotype_dic[barcode]] for barcode in umap_coordinates.index]
    ax.scatter(umap_coordinates['UMAP_1'].values, umap_coordinates['UMAP_2'].values, c=col, s=1)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    # ax.set_title('UMAP visualization of the macro cell types')
    for celltype in colors_dic: # print the macro types
        if celltype != 'nan':
            if celltype == 'T':
                pos = (-1, 3)
            else:
                pos = (0, 0)
                sub_cells = df_macro.loc[df_macro['Label'] == celltype].index
                bbox = macro_bbox[celltype]
                n_cells = 0
                for barcode in sub_cells:
                    cell_pos = (umap_coordinates['UMAP_1'][barcode], umap_coordinates['UMAP_2'][barcode])
                    if cell_pos[0] >= bbox[0][0] and cell_pos[0] <= bbox[0][1]:
                        if cell_pos[1] >= bbox[1][0] and cell_pos[1] <= bbox[1][1]:
                            pos = (pos[0] + cell_pos[0], pos[1] + cell_pos[1])
                            n_cells += 1
                pos = (pos[0] / n_cells, pos[1] / n_cells)

            if celltype == 'B':
                pos = (pos[0], pos[1]-0.5)

            ax.text(*pos, celltype, horizontalalignment='center', verticalalignment='center')

    legend_elements = [ Patch(facecolor=colors_dic[celltype], label=celltype if celltype != 'nan' else 'NA') for celltype in colors_dic]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    ax.set_aspect((ax.get_xlim()[1] - ax.get_xlim()[0])/(ax.get_ylim()[1] - ax.get_ylim()[0]))

    fig.tight_layout(h_pad=1)

    plt.show()

    return

#-------------------------------------------------------------------------------
def plot_LEF1():

    # load the UMAP coordinates
    file_name = '../../dataset/Imagine/umap_coordinates.csv'
    umap_coordinates = pd.read_csv(file_name, index_col = 0)

    # import the cell types
    file_name = '../../dataset/Imagine/cell_types.csv'
    df_cell_types = pd.read_csv(file_name, index_col=0)
    df_cell_types.rename(columns={'cellType_final': 'Label'}, inplace=True)
    print(df_cell_types)

    # load the normalized matrix
    file_name = '../../dataset/Imagine/normalized_matrix.csv'
    df = pd.read_csv(file_name, index_col = 0)
    df = df.T


    # plot the LEF1 expression
    cell_indexes_sorted = list(umap_coordinates.index.values)
    LEF1_values = df.loc[cell_indexes_sorted]['LEF1'].values
    LEF1_values_dic = {cell_indexes_sorted[ind]:LEF1_values[ind] for ind in range(len(cell_indexes_sorted))}


    cell_indexes_sorted.sort(key=lambda ind:LEF1_values_dic[ind])

    LEF1_values = df.loc[cell_indexes_sorted]['LEF1'].values

    fig, ax = plt.subplots()
    umap_coordinates_bis = umap_coordinates.loc[cell_indexes_sorted]
    ax.scatter(umap_coordinates_bis['UMAP_1'], umap_coordinates_bis['UMAP_2'], marker='o', s=4, c=LEF1_values)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('LEF 1')
    ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))


    y_dic = {ind:np.random.normal(0,1) for ind in df.index}

    fig, axs = plt.subplots(1,3)

    ax = axs[0]
    ax.scatter(np.random.normal(0,1, len(LEF1_values)), LEF1_values, marker='x')

    CD8_naives = [ind for ind in df_cell_types.index if df_cell_types['Label'][ind] == 'CD8-naive']
    CD4_naives = [ind for ind in df_cell_types.index if df_cell_types['Label'][ind] == 'CD4-naive']

    CD4_CD8_immatures = CD8_naives + CD4_naives

    CD8 = [ind for ind in df_cell_types.index if df_cell_types['Label'][ind] == 'CD8']
    CD4 = [ind for ind in df_cell_types.index if df_cell_types['Label'][ind] == 'CD4']

    CD4_CD8 = CD4+CD8

    others = [ind for ind in df.index if not ind in CD4_CD8 and not ind in CD4_CD8_immatures]

    ax = axs[1]
    # violin plot for each gene
    # parts = ax.violinplot(df.loc[CD4_CD8_immatures]['LEF1'].values, showmeans=False, showmedians=True)
    parts = ax.violinplot([df.loc[CD4_CD8]['LEF1'].values, df.loc[CD4_CD8_immatures]['LEF1'].values, df.loc[others]['LEF1']], showmeans=False, showmedians=True)

    for pc in parts['bodies']:
        pc.set_facecolor('green')
        pc.set_edgecolor('black')
        pc.set_alpha(0.5)

    ax.scatter([y_dic[index] for index in CD4_CD8_immatures], [df['LEF1'][index] for index in CD4_CD8_immatures], marker='x', c='red')
    ax.set_ylim(axs[0].get_ylim())
    ax.set_xlim(axs[0].get_xlim())

    ax = axs[2]
    ax.scatter([y_dic[index] for index in CD4_CD8], [df['LEF1'][index] for index in CD4_CD8], marker='x', c='green')
    ax.set_xlim(axs[0].get_xlim())
    ax.set_ylim(axs[0].get_ylim())

    plot_cell_types()

    return

#-------------------------------------------------------------------------------
def umap_visualisation():

    plot_cell_types()

    plot_cell_macrotypes()

# umap_visualisation()

plot_LEF1()
