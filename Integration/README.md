# Notes on .mex files

My Intel Fortran Compiler (*Intel Parallel Studio XE Cluster Edition for Windows*) is expired. So that I could not re-build *.mex functions from Fortran codes in this folder, I just copied them from the previous version. However, they work fine for now.

We do not have the mex version of the NearlySingularIntegrals_KxpN function. We will use MATLAB version for now; unfortunately it will slow down the computation a little bit.

This way of computing the integrals equations is ridiculous. I hope that we can find a better way to compute these element integrals efficiently and accurately.

- Note that using multiple dispatch for the element integrals in Julia language could be a way to deal with this mess.
