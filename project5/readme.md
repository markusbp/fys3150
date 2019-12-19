# Project 5, Code Content:

quantum_dots.cpp/hpp:
- Class implementation of 2 electron quantum dot system, in Harmonic Oscillator
potential, with different trial wavefunctions, for use with variational
Monte Carlo simulation.

psi1.cpp:
- Script for performing grid search for psi1 under variational parameter alpha

psi2.cpp:
- Script for running psi2 with desired variational parameters, alpha, beta

virial.cpp:
- Script to compute and save kinetic and potential energy for optimal psi2 wave func. Used to generate data for testing virial theorem for different HO frequencies

unit_tests.cpp:
- Very simple unit test to check acceptance ratio of both wavefuncs for optimal parameters alpha and beta, for different HO frequencies.

searc_psi2.cpp:
- 2D grid search of alpha and beta variational parameters for psi2, using MPI for parallelizing.

plot_device.py:
- script for plotting.

Makefile:
- makefile :9
