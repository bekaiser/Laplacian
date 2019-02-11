# MPI parallel computation of the two-dimensional Laplacian operator
This repository contains fortran 90 tutorial scripts that show how MPI can be used to divide a computational domain into subarrays that can be processed on separate ranks and then recombined. The tutorials will be especially helpful to those looking to parallelize finite difference operations or to parallelize spectral methods. 

The included scripts do the following:
1) mpi_scatter_gather_1d.f90: scatters a vector, adds 1, and then gathers the vector again on the master rank.
2) mpi_scatter_gather_1d.f90: scatters a 2d array, adds 1, and then gathers the 2d array again on the master rank.
3) mpi_cfd2.f90: computes the Laplacian of a signal (with an analytical solution) using second-order accurate finite differencing and outputs ASCII files for plotting. 
