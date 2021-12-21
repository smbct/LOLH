# LOLH: Learning Optimized Logical Hypothesis

## Background

LOLH is an Inductive Logic Programming framework based on Learning From Interpretation Transition (LFIT: https://github.com/Tony-sama/pylfit.git). Its goal is to learn logic rules from single cell gene expression data, on order to infer relations between the genes. LOLH performs a selection of genes based on a discrete optimization problem formulation. It is implemented in c++ and python 3.8.

## Dependencies

The following python packages are required to execute the code: scipy, numpy, pandas, networkx, python-louvain, umap-learn, matplotlib. pylfit is also required to perform a comparison between PRIDE and LOLH algorithms. The c++ code relies on openmp. The R script 'DE_testing.R' for differential expression (DE) analysis relies on the R package Seurat.

## Usage

The python source code is stored in the python folder and the c++ code is stored in the c++ folder. The c++ code needs to be compiled using the command 'make' from LOLH directory. The pre-processor definitions can be modified in Constants.hpp to activate parallel execution of the program (USE_OPENMP 1 and N_THREADS equal to the number of cores).

A complete demonstration of LOLH applied to a concrete dataset is proposed in the python script 'main_analysis.py'. The demonstration is intended to be executed on the dataset available at the following link: http://doi.org/10.5281/zenodo.4730807 . In order to be executed, the python script 'main_analysis.py' should be run from the folder 'python'. The files available in the dataset need to be placed in a 'dataset' folder, in LOLH parent directory. The parent directory should also contain the folders 'coexpression' and 'cell_classification', as illustrated below.

Part of the analysis associated to the selection of cells in the demo relies on the differential expression (DE) analysis tool in Seurat (https://satijalab.org/seurat/). The analysis can be performed using the R script 'DE_testing.R'.

## Examples

### artificial dataset:

Randomly generated LOLH instances.

### Imagine dataset:

A pbmc matrix obtained by Imagine institute (https://www.institutimagine.org/fr)

-> http://doi.org/10.5281/zenodo.4730807

## Acknowledgment:

## References:

Samuel Buchet, Francesco Carbone, Morgan Magnin, Mickaël Ménager, and Olivier Roux. 2021. Inference of Gene Networks from Single Cell Data through Quantified Inductive Logic Programming. In The 12th International Conference on Computational Systems-Biology and Bioinformatics (CSBio2021). Association for Computing Machinery, New York, NY, USA, 48–63. DOI:https://doi.org/10.1145/3486713.3486746

Ribeiro T., Folschette M., Magnin M., Roux O., Inoue K. (2018) Learning Dynamics with Synchronous, Asynchronous and General Semantics. In: Riguzzi F., Bellodi E., Zese R. (eds) Inductive Logic Programming. ILP 2018. Lecture Notes in Computer Science, vol 11105. Springer, Cham. https://doi.org/10.1007/978-3-319-99960-9_8

Tony Ribeiro, Maxime Folschette, Morgan Magnin, Katsumi Inoue. Polynomial Algorithm For Learning From Interpretation Transition. 1st International Joint Conference on Learning & Reasoning, Oct 2021, (virtual), Greece. ⟨hal-03347026⟩
