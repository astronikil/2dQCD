#! /bin/zsh



gfortran -O3 -W -o correl main.F checktolinear.F paulinit.F propagator.F c2pt.F ../addressing/addrc.F ../dslash/Doverlap.F ../dslash/Dwilson.F ../dslash/ferm.F ../invertor/invert_overlap.F ../invertor/invert_wilson_shifted.F ../gauge/wloops.F ../gauge/sun.F ../gauge/dummysmear.F ../gauge/unitar.F ../gauge/expsun.F ../util/io.F ../util/mydot.F ../util/gauss.F ../gauge/gtransform.F ../util/hconj.F -llapack -I ../headers
  




