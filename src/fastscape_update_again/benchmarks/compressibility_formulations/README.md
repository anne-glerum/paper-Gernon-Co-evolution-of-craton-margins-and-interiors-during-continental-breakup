This directory contains the benchmark cases for different compressibility
approximations discussed in:

Gassmöller, R., Dannberg, J., Bangerth, W., Heister, T., & Myhill, R. (2020). On formulations of compressible mantle convection. Geophysical Journal International, 221(2), 1264-1280.

All benchmarks require the additional ASPECT plugins that lie in the folder
'plugins' and need to be compiled into a shared library before running the
benchmark (see the manual here:
http://www.math.clemson.edu/~heister/manual.pdf#sec%3Abenchmark-run). The
individual benchmarks live in separate directories that each contain a shell
script 'run.sh' which runs the benchmark for different approximations and
resolutions. The run script in the top level directory executes all benchmarks
in turn. The 'figures' directory contains scripts for gnuplot that will
organize the results into convergence plots (no need to copy files around). The
resulting plots can be directly compared against the figures in the
publication.
