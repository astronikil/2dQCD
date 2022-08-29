SHELL = /bin/sh

BIN = correl

FC = gfortran
LDR = gfortran

FFLAGS  = -O3 -freal-4-real-8
LDFLAGS = -O3 -freal-4-real-8

.SUFFIXES:  .f .c .o

F_OBJ =  addrc.o ape.o boost.o c2pt.o expsun.o ferm.o invertor.o io.o smooth.o sun.o util.o wloops.o pauliinit.o dirprod.o main.o sources.o seqsource.o propagator.o seqpropagator.o invertandcontract.o printout.o

compile: $(BIN)

$(BIN): $(F_OBJ)
	$(LDR) $(LDFLAGS) -o $(BIN) $(F_OBJ) -llapack

.o.o:
	$(FC) $(FFLAGS) -c $*.o -llapack

clean:
	rm -f $(F_OBJ) $(C_OBJ) $(BIN) core *.i *.lst
