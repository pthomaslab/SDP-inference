# sdp-inf
Example codes for our paper **Moment-based parameter inference with error guarantees for stochastic reaction networks**.

---------------------------------
The file Schlogl_example.m contains MATLAB codes for obtaining the moment equations and moment matrices constraints, setting up constrained sets with these constraints via YALMIP, and optimising over the set to get bounds on parameters. These correspond to Fig.2(c) in the paper. One can use the example data of moment intervals (in the file Example_data.mat) to directly try the optimisation and obtain bounds on parameters.

The file PT_example.m has the same components, but additionally how to integrate partially-observed data with different species observed. These correspond to Fig.4 in the paper.

