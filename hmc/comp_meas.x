#! /bin/bash



gfortran -O3 -W -o hmc \
../hmc/main_meas.F \
../hmc/lat.F \
../hmc/metro.F \
../hmc/mold_leap.F \
../hmc/bosefo.F \
../hmc/energy.F \
../hmc/evol.F \
../hmc/fermfo.F \
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
../ritz_rajamani/submain.f \
../ritz_rajamani/ddag.F \
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
-llapack \
-I ../headers \
-I ../ritz_rajamani


