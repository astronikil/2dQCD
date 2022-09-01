#! /bin/bash



gfortran -O3 -W -o correl ../c2pt/main.F ../c2pt/checktolinear.F ../c2pt/paulinit.F ../c2pt/propagator.F ../c2pt/c2pt.F ../addressing/addrc.F ../dslash/Doverlap.F ../dslash/Dwilson.F ../dslash/ferm.F ../invertor/invert_overlap.F ../invertor/invert_wilson_shifted.F ../gauge/wloops.F ../gauge/sun.F ../gauge/dummysmear.F ../gauge/unitar.F ../gauge/expsun.F ../util/io.F ../util/mydot.F ../util/gauss.F ../gauge/gtransform.F ../util/hconj.F -llapack -I ../headers
  




