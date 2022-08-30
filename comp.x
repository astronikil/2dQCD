#! /bin/zsh



gfortran -O3 -W -o su2 dslash/Doverlap.F dslash/Dwilson.F addressing/addrc.F hmc/bosefo.F gauge/dummysmear.F hmc/energy.F hmc/evol.F gauge/expsun.F dslash/ferm.F hmc/fermfo.F invertor/invert_overlap.F invertor/invert_wilson_shifted.F util/io.F hmc/lat.F hmc/main.F hmc/metro.F hmc/mold_leap.F gauge/sun.F gauge/unitar.F util/mydot.F rvector/phiinit.F rvector/gnoise.F rvector/ranv.F gauge/wloops.F gauge/sungenerators.F util/gauss.F -llapack -I ./headers



#module load openmpi/3.1.4/intel-2018

#mpifort -O3 -xCORE-AVX2 -fma -align array64byte -finline-functions  -o su2 util.F ferm.F lat.F main.F mold.F sun.F io.F addrc.F vector.F congrad_simple.F wloops.F metro.F evol.F  unitar.F -L../lapack-3.9.0/ -llapack


#gfortran -O3 -W -o su2 util.F ferm.F lat.F main.F mold.F sun.F io.F addrc.F vector.F invertor.F wloops.F metro.F evol.F  unitar.F  expsun.F -llapack

#gfortran -O3 -W -o su2 util.F ferm.F lat.F main.F mold.F sun.F io.F addrc.F vector.F invertor.F wloops.F metro.F evol.F  unitar.F  expsun.F dummysmear.F Doverlap.F Dwilson.F invert_overlap.F invert_wilson_shifted.F -llapack
