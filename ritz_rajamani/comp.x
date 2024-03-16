#! /bin/bash


#gfortran -O3 -o suneig ../su2J/fundtoJ.f ../su2J/su2alg.f ../ritz_rajamani/check.F ../ritz_rajamani/ddag.F ../ritz_rajamani/random.f ../ritz_rajamani/mathtools.f ../ritz_rajamani/sort.f ../ritz_rajamani/spectrum.f ../ritz_rajamani/startritz.f ../ritz_rajamani/ritz.f ../ritz_rajamani/diagritz.f ../ritz_rajamani/GS.f ../dslash/ferm.F ../dslash/Dwilson.F ../dslash/Doverlap.F ../invertor/invert_overlap_shifted.F ../invertor/invert_wilson_shifted.F ../addressing/addrc.F ../util/io.F ../gauge/wloops.F  ../gauge/dummysmear.F ../ritz_rajamani/miscroutines.f ../ritz_rajamani/ipr.f ../util/mydot.F ../gauge/sun.F ../gauge/unitar.F ../gauge/sungenerators.F ../gauge/expsun.F ../util/gauss.F ../gauge/wloops_smear.F -I ../headers -I ../../su2J -llapack

#gfortran -O3 -o suneig ../su2J/fundtoJ.f ../su2J/su2alg.f ../ritz_rajamani/check.F ../ritz_rajamani/ddag.F ../ritz_rajamani/random.f ../ritz_rajamani/mathtools.f ../ritz_rajamani/sort.f ../ritz_rajamani/spectrum.f ../ritz_rajamani/startritz.f ../ritz_rajamani/ritz.f ../ritz_rajamani/diagritz.f ../ritz_rajamani/GS.f ../dslash/ferm.F ../dslash/Dwilson.F ../dslash/Doverlap.F ../invertor/invert_overlap_shifted.F ../invertor/invert_wilson_shifted.F ../addressing/addrc.F ../util/io.F ../gauge/wloops.F ../gauge/wloops_smear.F ../gauge/stout.F ../gauge/expsun_stout.F ../ritz_rajamani/miscroutines.f ../ritz_rajamani/ipr.f ../util/mydot.F ../gauge/sun.F ../gauge/unitar.F ../gauge/sungenerators.F ../gauge/expsun.F ../util/gauss.F -I ../headers -I ../../su2J -llapack

#compile for zero mode in z2 flux background 

gfortran -O3 -o suneig ../su2J/fundtoJ.f ../su2J/su2alg.f ../ritz_rajamani/check.F ../ritz_rajamani/ddag.F ../ritz_rajamani/random.f ../ritz_rajamani/mathtools.f ../ritz_rajamani/sort.f ../ritz_rajamani/spectrum.f ../ritz_rajamani/startritz.f ../ritz_rajamani/ritz.f ../ritz_rajamani/diagritz.f ../ritz_rajamani/GS.f ../dslash/ferm.F ../dslash/Dwilson.F ../dslash/Doverlap.F ../invertor/invert_overlap_shifted.F ../invertor/invert_wilson_shifted.F ../addressing/addrc.F ../util/io.F ../gauge/wloops.F ../gauge/wloops_smear.F ../gauge/stout.F ../gauge/expsun_stout.F ../ritz_rajamani/miscroutines.f ../ritz_rajamani/ipr.f ../util/mydot.F ../gauge/sun.F ../gauge/unitar.F ../gauge/sungenerators.F ../gauge/expsun.F ../util/gauss.F ../util/readinit.F -I ../headers  -llapack


