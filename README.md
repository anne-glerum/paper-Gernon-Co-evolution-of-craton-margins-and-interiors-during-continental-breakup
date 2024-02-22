This repository provides the input file, ASPECT installation details and additional ASPECT plugins used for the paper

*Co-evolution of craton margins and interiors during continental breakup*

by

Gernon, T. M.,
Hincks, T. K.,
Brune, S.,
Braun, J.
Jones, S. M.,
Keir, D.,
Cunningham, A. and
Glerum, A.

currently under review in Nature.

# Documentation
The numerical simulations presented in this paper were run with the geodynamics code ASPECT (https://aspect.geodynamics.org/).

## ASPECT version
The ASPECT input files provided in this repository correspond to commit a1f0aa5 of the ASPECT branch 

https://github.com/EstherHeck/aspect/tree/fastscape_update_again

This branch is built on commit 84d40e7 of the ASPECT main version,
which can be found at https://github.com/geodynamics/aspect. A copy of a1f0aa5 can be found in the
folder /src.

## Additional ASPECT plugins
For the initial model conditions, we used the ASPECT plugins in the folder /plugins. 
The file CMakeLists.txt can be used to install these plugins as shared libraries
against your ASPECT installation.

## ASPECT input file
The ASPECT input files can be found in the folder /prms. The file names indicate what simulation in the paper
they correspond to.

## Installation details
ASPECT was built using the underlying library deal.II 10.0.0-pre (master, 18a3861)
on the German HLRN cluster Lise. deal.II used:
* 32 bit indices and vectorization level 3 (512 bits)
* Trilinos 12.18.1
* p4est 2.2.0
