#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.patches import Patch

def UMAP_visualisation():

    # read the normalized matrix (logcpm)
    filename = '../../dataset/Imagine/normalized_dataset.csv'
    df = pd.read_csv(filename, index_col=0)
    df = df.T

    # load the UMAP cells coordinates
    df_coordinates = pd.read_csv('../../dataset/UMAP_coordinates.csv', index_col=0)

    # read the seurat cell types
    df_celltypes = pd.read_csv('../../dataset/cell_types.csv', index_col = 0)
    df_celltypes.rename(columns={'cellType_final': 'Label'}, inplace=True)
    celltype_dic = {str(barcode):str(df_celltypes['Label'][barcode]) for barcode in df_celltypes.index}
    print(df_celltypes)

    # read the seurat "macro" types
    df_macro = pd.read_csv('../../dataset/cell_types_macro.csv', index_col = 0)
    df_macro.rename(columns={'cellType_macro': 'Label'}, inplace=True)
    macro_types = [str(label) for label in np.unique(list(df_macro['Label'].values))]
    macrotype_dic = {str(barcode):str(df_macro['Label'][barcode]) for barcode in df_macro.index}
    print(df_macro)



    ###############################################################################
    # plot the UMAP with the cell types
    ###############################################################################

    colors = [(248, 118, 109), (232, 133, 38), (211, 146, 0), (183, 159, 0), (147, 170, 0), (94, 179, 0)]
    colors += [(0, 186, 56), (0, 191, 116), (0, 193, 159), (0, 191, 196), (0, 185, 227), (0, 173, 250)]
    colors += [(97, 156, 255), (174, 135, 255), (219, 114, 251), (245, 100, 227), (255, 97, 195), (255, 105, 156)]
    for index in range(len(colors)): # thanks matplotlib..
        col = colors[index]
        colors[index] = (col[0]/255., col[1]/255., col[2]/255.)

    cell_types = ['CD4-naive', 'CD14', 'CD4', 'NK', 'CD8', 'CD16', 'B-naive', 'Treg']
    cell_types += ['gdT', 'CD8-naive', 'B', 'T_undef', 'cDC', 'platelet_cont', 'CD14-PPBP', 'HSC', 'pDC', 'undefined']

    colors_dic = {cell_types[index]:colors[index] for index in range(len(cell_types))}

    fig, ax = plt.subplots()
    col = [colors_dic[celltype_dic[barcode]] for barcode in df_coordinates.index]
    ax.scatter(df_coordinates['UMAP_1'].values, df_coordinates['UMAP_2'].values, c=col, s=1)

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('UMAP visualization of the cell types')

    # print the cell types
    for celltype in cell_types:

        if celltype != 'undefined':

            if celltype == 'T_undef':
                pos = (1.15, -0.15)
            elif celltype == 'Treg':
                pos = (-0.3, 1.3)
            else:
                pos = (0, 0)
                sub_cells = df_celltypes.loc[df_celltypes['Label'] == celltype].index
                for barcode in sub_cells:
                    pos = (pos[0] + df_coordinates['UMAP_1'][barcode], pos[1]+df_coordinates['UMAP_2'][barcode])
                pos = (pos[0] / len(sub_cells), pos[1] / len(sub_cells))
            ax.text(*pos, celltype, horizontalalignment='center', verticalalignment='center')

    legend_elements = [ Patch(facecolor=colors_dic[celltype], label=celltype) for celltype in colors_dic]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    fig.tight_layout(h_pad=1)



    ###############################################################################
    # plot the UMAP with the cell types
    ###############################################################################

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
    col = [colors_dic[macrotype_dic[barcode]] for barcode in df_coordinates.index]
    ax.scatter(df_coordinates['UMAP_1'].values, df_coordinates['UMAP_2'].values, c=col, s=1)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('UMAP visualization of the macro cell types')
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
                    cell_pos = (df_coordinates['UMAP_1'][barcode], df_coordinates['UMAP_2'][barcode])
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

    fig.tight_layout(h_pad=1)

    plt.show()

# UMAP_visualisation()
