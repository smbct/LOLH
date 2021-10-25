#!/usr/bin/python3

import numpy as np
import pandas as pd

#-------------------------------------------------------------------------------
def print_lines(df):

    for index in df.index:

        # print(df.loc[index])

        line = index + ' & '

        ind = 0
        for val in df.loc[index]:
            line += '$'
            # line += str(val)
            val_str = np.format_float_scientific(val, precision = 2, exp_digits=3)

            if val_str[-3:] == '000':
                val_str = val_str[:-5]

            line += val_str

            line += '$'
            if ind < df.shape[1]-1:
                line += ' & '
            ind += 1
        line += ' \\\\ \hline'
        print(line)

    return

# import a csv file from DE testing and output a tex formatted text

# filename = 'NK_classification_markers.csv'
# df = pd.read_csv(filename, index_col=0)
# df = df.loc[df.index[:30]]

filename = 'c6_markers.csv'

# filename = 'c5_DE.csv'

df = pd.read_csv(filename, index_col=0)
df_sub = df.loc[df.index[:20]]
print_lines(df_sub)
print('\n\n')
df_sub = df.loc[df.index[-20:]]
print_lines(df_sub)



# print(df.head(), '\n\n')
