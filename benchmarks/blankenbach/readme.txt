Blankenback benchmark

This folder allows running Blankenbach case 1a, 2b, 2a, and 2b from the paper

Blankenbach, B., et al. "A benchmark comparison for mantle convection codes."
Geophysical Journal International 98.1 (1989): 23-38.

The reference values (given in caseXX_reference.stat) are from that paper (see
Table 9). When you run the benchmarks by typing in "make" after building the
plugin in the plugin/ subfolder, you will be eventually presented with the
output below. Pleas note that these computations take a long time to reach
steady state, especially on finer meshes.

See the comment in the "Makefile" in this folder for more options for running
the computations.

Output:
# Nu           Vrms           name:
4.87452472e+00 4.28499932e+01 case1a_ref3.stat
4.88459738e+00 4.28656773e+01 case1a_ref4.stat
4.88442981e+00 4.28650353e+01 case1a_ref5.stat
4.88441081e+00 4.28649624e+01 case1a_ref6.stat
4.88440900e+00 4.28649470e+01 case1a_reference.stat
1.03258810e+01 1.95279365e+02 case1b_ref3.stat
1.04949361e+01 1.93087841e+02 case1b_ref4.stat
1.05339127e+01 1.93214560e+02 case1b_ref5.stat
1.05339404e+01 1.93214791e+02 case1b_ref6.stat
1.05340950e+01 1.93214540e+02 case1b_reference.stat
2.24226685e+01 8.89265354e+02 case1c_ref3.stat
2.12315842e+01 8.39487685e+02 case1c_ref4.stat
2.18584328e+01 8.33321113e+02 case1c_ref5.stat
2.19710401e+01 8.33972142e+02 case1c_ref6.stat
2.19724650e+01 8.33989770e+02 case1c_reference.stat
7.27920586e+00 4.39260072e+02 case2a_ref3.stat
1.02189887e+01 4.64971807e+02 case2a_ref4.stat
1.01930353e+01 4.75819129e+02 case2a_ref5.stat
1.00711673e+01 4.80164721e+02 case2a_ref6.stat
1.00660000e+01 4.80433400e+02 case2a_reference.stat
6.44510878e+00 1.68574818e+02 case2b_ref3.stat
6.94481365e+00 1.69976109e+02 case2b_ref4.stat
6.92988247e+00 1.71554254e+02 case2b_ref5.stat
6.92964921e+00 1.71744538e+02 case2b_ref6.stat
6.92990000e+00 1.71755000e+02 case2b_reference.stat

steps needed:
72 output-case1a_ref3/statistics
128 output-case1a_ref4/statistics
232 output-case1a_ref5/statistics
470 output-case1a_ref6/statistics
210 output-case1b_ref3/statistics
347 output-case1b_ref4/statistics
690 output-case1b_ref5/statistics
1368 output-case1b_ref6/statistics
830 output-case1c_ref3/statistics
1916 output-case1c_ref4/statistics
3397 output-case1c_ref5/statistics
6728 output-case1c_ref6/statistics
1857 output-case2a_ref3/statistics
3088 output-case2a_ref4/statistics
6638 output-case2a_ref5/statistics
13593 output-case2a_ref6/statistics
4076 output-case2b_ref3/statistics
8404 output-case2b_ref4/statistics
17034 output-case2b_ref5/statistics
34164 output-case2b_ref6/statistics
