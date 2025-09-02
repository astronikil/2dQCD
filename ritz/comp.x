#! /bin/bash


gfortran -O3 -o suneig ../ritz/main.f ../ritz/ddag.F ../ritz/random.f ../ritz/mathtools.f ../ritz/sort.f ../ritz/spectrum.f ../ritz/startritz.f ../ritz/ritz.f ../ritz/diagritz.f ../ritz/GS.f ../dslash/ferm.F ../dslash/Dwilson.F ../dslash/Doverlap.F ../invertor/invert_wilson_shifted.F ../addressing/addrc.F ../util/io.F ../gauge/wloops.F  ../gauge/dummysmear.F ../ritz/miscroutines.f ../ritz/ipr.f ../util/mydot.F ../gauge/sun.F ../gauge/unitar.F ../gauge/sungenerators.F ../gauge/expsun.F ../util/gauss.F -I ../headers -llapack


