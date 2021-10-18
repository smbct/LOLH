#!/bin/bash

# download the single cell pbmc matrix from Imagine Institute

# zenodo link to the files
link='https://zenodo.org/record/4730807/files/'


celltypes_file='cell_types.csv'
macrotypes_file='cell_types_macro.csv?download=1'
umap_file='UMAP_coordinates.csv'
raw_matrix_file='raw_dataset.csv'
normalized_matrix_file='normalized_dataset.csv'

cd '../../dataset/' && mkdir 'Imagine'

cd 'Imagine'

mkdir 'coexpression'
mkdir 'cell_classification'

# download all the files from zenodo

wget $link$celltypes_file'?download=1' -O 'cell_types.csv'
wget $link$macrotypes_file'?download=1' -O 'cell_types_macro.csv'
wget $link$umap_file'?download=1' -O 'umap_coordinates.csv'
wget $link$normalized_matrix_file'?download=1' -O 'raw_matrix.csv'
wget $link$normalized_matrix_file'?download=1' -O 'normalized_matrix.csv'


# discretization of the normalized matrix
cd '../../example/Imagine/' && python3 'discretization.py'
