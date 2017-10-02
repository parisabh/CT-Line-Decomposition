#README

The source codes are developed by Parisa Babaheidarian. 

The provided functions decompose sinogram measurements in CT imaging for each path into a given basis subspace. The number of energy bin measurements considered for each path is set to 10 as default value. The codes are developed for three choices of basis function including PCA, photo-electric and compton, and SPECK basis function. We used PCA with 5 dimension, PEC with 2, and SPECK with 10 dimensions.

The Jac* files compute the Jacobian matrices for the decomposition cost function and for the given basis function. 

The *LineDecomposition* source files include the MEX functions as for the MATLAB wrappers.

The driver functions call the L-BFGS-B optimization routine. These routines can be downloaded online at https://github.com/stephenbeckr/L-BFGS-B-C.

To run the MEX file in MATLAB, create a folder /src and add the source and header files to this folder the run the following command:

 mex LineDecomposition.cpp lbfgsb.c linesearch.c subalgorithms.c print.c linpack.c miniCBLAS.c timer.c

The source functions for other basis transforms can be run in a similar way.

The source code gLineDecompositionSparse.cpp adds the L1 sparsity constraint into the cost function.