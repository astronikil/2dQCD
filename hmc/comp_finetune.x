#! /bin/bash

module load gcc/9.2.0
module load cpu/0.15.4
module load netlib-lapack/3.8.0


gfortran -O3 -o hmc \
../hmc/main_finetune.F \
../hmc/lat.F \
../hmc/metro.F \
../hmc/adapt_tune.F \
../hmc/mold_leap.F \
../hmc/bosefo.F \
../hmc/energy.F \
../hmc/evol.F \
../hmc/fermfo.F \
../hmc/dbydaniso.F \
../dslash/Doverlap.F \
../dslash/Dwilson.F \
../dslash/ferm.F \
../addressing/addrc.F \
../gauge/dummysmear.F \
../gauge/expsun.F \
../gauge/sun.F \
../gauge/unitar.F \
../gauge/wloops.F \
../gauge/sungenerators.F \
../invertor/invert_overlap.F \
../invertor/invert_wilson_shifted.F \
../util/io.F \
../util/mydot.F \
../util/gauss.F \
../rvector/phiinit.F \
../rvector/gnoise.F \
../rvector/ranv.F \
../rvector/z2noise.F \
../ritz_rajamani/submain.F \
../ritz_rajamani/DDAG.F \
../ritz_rajamani/random.f \
../ritz_rajamani/mathtools.f \
../ritz_rajamani/sort.f \
../ritz_rajamani/GS.f \
../ritz_rajamani/spectrum.f \
../ritz_rajamani/startritz.f \
../ritz_rajamani/ritz.f \
../ritz_rajamani/diagritz.f \
../ritz_rajamani/miscroutines.f \
../ritz_rajamani/ipr.f \
../c2pt/submain.F \
../c2pt/checktolinear.F \
../c2pt/paulinit.F \
../c2pt/propagator.F \
../c2pt/c2pt.F \
-llapack \
-I ../headers \
-I ../ritz_rajamani


