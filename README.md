# 2dQCD

This is a fortran based code suite for SU(Nc) gauge theory coupled to fermions in two-dimensions.

The code suite consists of subroutines that can be put together to construct 
1. Hybrid Monte Carlo 
2. Dirac eigenvalue measurement
3. Meson spectroscopy -- 2-point, 3-point function

The header files are contained in the directory ./headers/ . The file size.h contains all the parameters.

#Gauge action

The current implementation uses an unimproved Wilson gauge action
The coupling is the dimensionless 2d area, \lambda, that is related to lattice coupling \beta = L^2/\lambda.

#Fermions

The current implementation has 
1. Wilson-Dirac
2. Overlap
3. The hermitian operators of the above.

#Invertors
The inversion is via conjugate gradient. The current version has
1. single inversion of Wilson-Dirac, D^\dagger D
2. Single inversion of hermitian-overlap^2
3. Multimass inversion of hermitian-Dirac^2

The stopping criterion & max steps for for Wilson-Dirac are specified in block 
common /param_sqrt/ delitsq, nitsq
The stopping criterion & max steps for for Overlap are specified in block 
common /param_overlap/ delitov, nitov
