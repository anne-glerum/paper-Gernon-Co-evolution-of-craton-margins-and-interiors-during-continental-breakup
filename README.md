# paper-Gernon-Kimberlite-ascent-by-rift-driven-disruption

This repository provides the input file, ASPECT installation details and additional ASPECT plugins used for the paper

*Rift-induced disruption of cratonic keels drives kimberlite volcanism*

by

Gernon, T. M. 
Jones, S. M.
Brune, S.
Hincks, T. K.
Palmer, M. R.
Schumacher, J. C.
Primiceri, R. M.
Field, M.
Griffin, W. L. 
O'Reilly, S. Y.
Keir, D.
Spencer, C. J.
Merdith, A. S.
Glerum, A.

in Nature (https://doi.org/10.1038/s41586-023-06193-3).

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
The ASPECT input file can be found in the folder /prms.

## Installation details
ASPECT was built using the underlying library deal.II 10.0.0-pre (master, 18a3861)
on the German HLRN cluster Lise.
