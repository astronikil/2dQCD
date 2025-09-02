#! /bin/bash


gfortran -O3 -o hmc \
../hmc/main_meas2.F \
../hmc/lat.F \
../hmc/metro.F \
../hmc/adapt.F \
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
../ritz/submain.f \
../ritz/ddag.F \
../ritz/random.f \
../ritz/mathtools.f \
../ritz/sort.f \
../ritz/GS.f \
../ritz/spectrum.f \
../ritz/startritz.f \
../ritz/ritz.f \
../ritz/diagritz.f \
../ritz/miscroutines.f \
../ritz/ipr.f \
../c2pt/submain.F \
../c2pt/checktolinear.F \
../c2pt/paulinit.F \
../c2pt/propagator.F \
../c2pt/c2pt.F \
-llapack \
-I ../headers \
-I ../ritz


