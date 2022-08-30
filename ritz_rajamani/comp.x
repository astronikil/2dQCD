#! /bin/bash


gfortran main.f ddag.F random.f mathtools.f sort.f spectrum.f startritz.f ritz.f diagritz.f GS.f ../dslash/ferm.F ../dslash/Dwilson.F ../dslash/Doverlap.F ../invertor/invert_wilson_shifted.F ../addressing/addrc.F ../util/io.F ../gauge/wloops.F  ../gauge/dummysmear.F miscroutines.f ipr.f ../util/mydot.F ../gauge/sun.F ../gauge/unitar.F ../gauge/sungenerators.F ../gauge/expsun.F ../util/gauss.F -I ../headers -llapack


