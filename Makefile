SHELL = /bin/sh

FC = gfortran

FFLAGS  = -O3 #-fdefault-real-8 #-fopenmp

.SUFFIXES:  .F .f .c .o


compile: $(BIN)

slash/Doverlap.F dslash/Dwilson.F addressing/addrc.F hmc/bosefo.F gauge/dummysmear.F hmc/energy.F hmc/evol.F gauge/expsun.F dslash/ferm.F hmc/fermfo.F invertor/invert_overlap.F invertor/invert_wilson_shifted.F util/io.F hmc/lat.F hmc/main.F hmc/metro.F hmc/mold_leap.F gauge/sun.F gauge/unitar.F util/mydot.F rvector/phiinit.F rvector/gnoise.F rvector/ranv.F gauge/wloops.F gauge/sungenerators.F util/gauss.F

hmc: headers/size.h  hmc/main.o
	$(FC) $(FFLAGS) -o hmc   -llapack

.F.o:
	$(FC) $(FFLAGS) -c $*.F -llapack

clean:
	rm -f $(F_OBJ) $(C_OBJ) $(BIN) core *.i *.lst *.optrpt
