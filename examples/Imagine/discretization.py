#!/usr/bin/python

import numpy as np
import scipy.stats
import pandas as pd

#------------------------------------------------------------------------------
def discretization(filename, filename_discrete):

    df = pd.read_csv(filename, index_col=0)
    df = df.T

    # print(df)
    # print(df.shape)


    # dictionary gene to index
    ind = 0
    gene_ind = {}
    for gene in df.columns:
        gene_ind[gene] = ind
        ind += 1

    # creation of the discrete dataframe

    binary_data = [ [ 0 for col in range(df.shape[1])] for row in range(df.shape[0])]

    # computation of the gene medians
    medians = [ np.median(df[elt].values) for elt in df.columns ]

    # based on the semi-manual disctretization (see pdf)

    col_ind = 0
    for gene in df.columns:

        if col_ind%1000 == 0:
            print(col_ind, ' genes done over ', df.shape[1])


        median = medians[gene_ind[gene]]

        ############################################################
        if median <= 0.1: # one threshold

            # predefined manual thresholds
            predef = {'CD79A':1.0, 'IGHM':1.6, 'IGKC':2.0, 'VCAN':1.2, 'GZMH':1.25, 'FCN1':1.2, 'AIF1':1.4, 'PSAP':1.7, 'CS73':1.3, 'FOS':1.7 , 'S100A8':1.6, 'CST7':1.0, 'GZMA':1.0, 'TYROBP':1.2, 'GNLY':1.5, 'S100A9':1.5, 'LYZ':2.1, 'CCL5':1.5, 'NKG7':1.9}

            skew = scipy.stats.skew(df[gene].values)
            var = np.var(df[gene].values)

            threshold = None

            # look for the unique threshold
            if skew >= 4:
                threshold = 0.1
            elif skew <= 4:
                if var <= 0.3:
                    if skew <= 2:
                        threshold = 0.4
                    else:
                        threshold = 0.9
                else: # var >= 0.3
                    if gene in predef:
                        threshold = predef[gene]
                    else:
                        threshold=0.9

            # discretization in all the cels
            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if value >= threshold:
                    binary_data[row_ind][col_ind] = 1
                else:
                    binary_data[row_ind][col_ind] = 0
                row_ind += 1



        ############################################################
        elif median >= 0.6 and median <= 0.8:

            skew = scipy.stats.skew(df[gene].values)
            var = np.var(df[gene].values)

            predef = {'SAT1':2.1, 'COTL1':2.0, 'IL7R':1.0, 'CTSS':2.2}

            threshold = None

            # look for the unique threshold
            if var <= 0.45:
                if skew >= 0.4:
                    threshold = 0.9
                else:
                    threshold = 1.2
            else:
                if gene in predef:
                    threshold = predef[gene]
                else:
                    threshold = 0.9

            # discretization in all the cels
            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if value >= threshold:
                    binary_data[row_ind][col_ind] = 1
                else:
                    binary_data[row_ind][col_ind] = 0
                row_ind += 1



        ############################################################
        elif median >= 1.0 and median <= 1.2:

            # upper discretisation

            median = medians[gene_ind[gene]]

            if gene != 'LTB':
                threshold = 1.3*np.mean([val for val in df[gene].values if val > median+10**-4])

            # discretization in all the cels
            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if gene == 'LTB':
                    if value >= 3.2:
                        binary_data[row_ind][col_ind] = 2
                    elif value >= 1.2:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                else:
                    if value >= threshold:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                row_ind += 1




        ############################################################
        elif median >= 1.3 and median <= 1.5:

            predef = ['FYB1', 'TSC22D3', 'ETS1', 'HCST', 'SRGN']

            median = medians[gene_ind[gene]]

            threshold2 = 1.3 * np.mean([val for val in df[gene] if val > median+10**-4])

            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if gene in predef: # one threshold
                    if value >= threshold2:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                elif gene == 'CD3E':
                    if value >= threshold2:
                        binary_data[row_ind][col_ind] = 2
                    elif value >= 0.9:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                else:
                    if value >= threshold2:
                        binary_data[row_ind][col_ind] = 2
                    elif value >= 0.4:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                row_ind += 1



        ############################################################
        elif median >= 1.6 and median <= 1.7:

            predef = ['CYBA', 'ZFP36L2', 'JUNB', 'VIM']
            median = medians[gene_ind[gene]]

            if gene != 'IL32' and gene != 'S100A9':
                threshold2 = 1.2 * np.mean([val for val in df[gene] if val > median+10**-4])

            row_ind = 0
            for barcode in df.index:

                value = df.loc[barcode,gene]

                if gene == 'IL32':
                    if value >= 3.5:
                        binary_data[row_ind][col_ind] = 2
                    elif value >= 1.2:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                elif gene == 'S100A9':
                    if value >= 3.5:
                        binary_data[row_ind][col_ind] = 2
                    elif value >= 0.2:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                elif gene in predef:
                    if value >= threshold2:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                else:
                    if value >= threshold2:
                        binary_data[row_ind][col_ind] = 2
                    elif value >= 0.3:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0

                row_ind += 1



        ############################################################
        elif median >= 1.7 and median <= 1.9:

            median = medians[gene_ind[gene]]
            threshold1 = None
            threshold2 = 1.2 * np.mean([val for val in df[gene] if val > median+10**-4])

            if gene == 'S100A4':
                threshold1 = 1.5
            elif gene in ['NPM1', 'MT-ND5', 'KLF2', 'CALM1']:
                threshold1 = 0.2
            else:
                threshold1 = 0.8

            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if value >= threshold2:
                    binary_data[row_ind][col_ind] = 2
                elif value >= threshold1:
                    binary_data[row_ind][col_ind] = 1
                else:
                    binary_data[row_ind][col_ind] = 0
                row_ind += 1



        ############################################################
        elif median >= 1.9 and median <= 2.0:

            median = medians[gene_ind[gene]]
            threshold1 = 0.9
            threshold2 = 1.2 * np.mean([val for val in df[gene] if val > median+10**-4])

            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if value >= threshold2:
                    binary_data[row_ind][col_ind] = 2
                elif value >= threshold1:
                    binary_data[row_ind][col_ind] = 1
                else:
                    binary_data[row_ind][col_ind] = 0
                row_ind += 1



        ############################################################
        elif median >= 2.0 and median <= 2.15:

            median = medians[gene_ind[gene]]
            threshold2 = 1.2 * np.mean([val for val in df[gene] if val > median+10**-4])

            if gene == 'CD52' or gene == 'SH3BGRL3':
                threshold1 = 0.2
            else:
                threshold1 = 0.9

            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if value >= threshold2:
                    binary_data[row_ind][col_ind] = 2
                elif value >= threshold1:
                    binary_data[row_ind][col_ind] = 1
                else:
                    binary_data[row_ind][col_ind] = 0
                row_ind += 1



        ############################################################
        else:

            median = medians[gene_ind[gene]]

            if gene != 'FTL' and gene != 'FTH1':
                threshold1 = 0.8*np.mean([val for val in df[gene] if val < median+10**-4])
                threshold3 = 1.12*np.mean([val for val in df[gene] if val > median+10**-4])

            predef1 = ['RPS25', 'RPS5', 'RPS6', 'RPS12', 'EEF1B2', 'RPL35A', 'RPL5', 'RPL30', 'RPS16', 'RPL9', 'RPS8', 'RPL32', 'RPS13', 'RPL22', 'RPLP0', 'RPL38', 'RPL18', 'RPL10A', 'RPS20', 'RPL11', 'RPL19', 'RPL13', 'RPLP2', 'RPS3', 'RPL21', 'RPS14', 'RPS28', 'RPL18-A', 'TPT1', 'RPL34', 'RPL10', 'RPL37', 'RPL41']

            predef2 = {'ACTB':3.6, 'RPS21':2.98, 'RPS3A':3.62, 'RPS27A':3.35, 'HLA-A':3.02, 'RPS18':3.3, 'HLA-C':3.2, 'B2M':4.62, 'RPL36':2.8, 'RPS29':3.35, 'RPS27':3.86}

            row_ind = 0
            for barcode in df.index:
                value = df.loc[barcode,gene]
                if gene == 'FTL':
                    if value >= 3.5:
                        binary_data[row_ind][col_ind] = 1
                    else:
                        binary_data[row_ind][col_ind] = 0
                else:
                    if gene in predef1:
                        threshold2 = median
                        if value >= threshold3:
                            binary_data[row_ind][col_ind] = 3
                        elif value >= threshold2:
                            binary_data[row_ind][col_ind] = 2
                        elif value >= threshold1:
                            binary_data[row_ind][col_ind] = 1
                        else:
                            binary_data[row_ind][col_ind] = 0
                    elif gene in predef2:
                        threshold2 = predef2[gene]
                        if value >= threshold3:
                            binary_data[row_ind][col_ind] = 3
                        elif value >= threshold2:
                            binary_data[row_ind][col_ind] = 2
                        elif value >= threshold1:
                            binary_data[row_ind][col_ind] = 1
                        else:
                            binary_data[row_ind][col_ind] = 0
                    else:
                        if value >= threshold3:
                            binary_data[row_ind][col_ind] = 2
                        elif value >= threshold1:
                            binary_data[row_ind][col_ind] = 1
                        else:
                            binary_data[row_ind][col_ind] = 0
                row_ind += 1


        col_ind += 1


    df_discrete = pd.DataFrame(binary_data, index = df.index, columns = df.columns)

    print(df_discrete)

    df_discrete.to_csv(filename_discrete)

discretization('../../dataset/Imagine/normalized_matrix.csv', '../../dataset/Imagine/discrete_matrix.csv')
