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


#How to run

1. First set lattice size and Nc in headers/size.h
2. You can change Zolotorev parameters in headers/zolo.h
3. Set number of eigenvalue measurements in ritz_rajamani/global.f
4. Go to the directory EXAMPLE.  Run ./comp.x to compile hmc+measurment routines. 
5. Copy the executable hmc and script run.x to another directory where you want to start running.
6. open script run.x, and set the coupling \lambda, and also check for other parameters.
7. Run ./run.x.
8. It will run, and produce plaquette and Polyakov loop measurement in data000. The file outrun000 will
   contain general runtime history, acceptance rates, etc.  The eigenvalue measurements will be in EIG/.
   Meson measurements in MES/. If variable iboolwrite==1, then gauge fields will be written in GAUGE/.

