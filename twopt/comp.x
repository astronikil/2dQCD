#! /bin/tcsh

#source /usr/local/etc/templates/.cshrc.el7-xeon

#ifort -O3 -xCORE-AVX2 -fma -align array64byte -finline-functions -r8 -o correl addrc.F smooth2d.F boost.F c2pt.F expsun.F ferm.F invertor.F io.F smooth.F sun.F util.F wloops.F pauliinit.F dirprod.F main.F sources.F seqsource.F propagator.F seqpropagator.F invertandcontract.F printout.F  fpi.F -L../lapack-3.9.0/ -llapack

module load openmpi/3.1.4/intel-2018

mpifort -O3 -xCORE-AVX2 -fma -align array64byte -finline-functions -r8 -o correl addrc.F smooth2d.F boost.F c2pt.F expsun.F ferm.F invertor.F io.F smooth.F sun.F util.F wloops.F pauliinit.F dirprod.F main.F sources.F seqsource.F propagator.F seqpropagator.F invertandcontract.F printout.F  fpi.F -L../lapack-3.9.0/ -llapack
