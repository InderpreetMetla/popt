Pseudospectral OPTimiser (POPT)

This software was developed as part of my undergraduate thesis at The University of Queensland. POPT is distributed under the MIT License. 

This software is able to solve continous time, multiple-phase trajectory optimisation problems 
using Legendre-Gauss-Radau orthogonal collocation methods.

POPT uses IPOPT as its NLP solver. Pre-compiled mex files for MATLAB are required to be obtained from: 

https://www.coin-or.org/download/binary/Ipopt/

IPOPT version 3.11.8 is recommended. The IPOPT mex files are to be placed inside the popt/nlp/ipopt directory.

POPT is able to use automatic differentiation to carry out derivative calculations. The open-source program MatlabAutoDiff developed by martinResearch has been interfaced with the solver and is able from: https://github.com/martinResearch/MatlabAutoDiff. Download this program and add it to your MATLAB path if you wish to use automatic differentiation. A finite-difference technique and a complex-step differentiation technique are in-built in POPT and require no other third-party programs. 

A user's guide for the software is available in the popt/docs directory.
