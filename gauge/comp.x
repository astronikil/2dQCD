#! /bin/bash


gfortran -O3 -o smplaq ../gauge/check.F ../addressing/addrc.F ../util/io.F ../gauge/wloops.F ../gauge/wloops_smear.F ../gauge/stout.F ../util/mydot.F ../gauge/sun.F ../gauge/unitar.F ../gauge/sungenerators.F ../gauge/expsun.F ../gauge/expsun_stout.F ../util/gauss.F -I ../headers  -llapack



