#! /bin/bash



gfortran -O3 -W -o hmc \
../hmc/lat.F \
../hmc/main.F \
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
-llapack \
-I ../headers


