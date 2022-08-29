SHELL = /bin/sh

BIN = correl

FC = gfortran
LDR = gfortran

FFLAGS  = -O3 -fdefault-real-8 #-fopenmp
LDFLAGS = -O3 -fdefault-real-8 #-fopenmp

.SUFFIXES:  .f .c .o

F_OBJ =  utilu1.o ferm.o main.o u1.o io.o addrc.o invertor.o wloops.o hyp.o clover.o c2pt.o

compile: $(BIN)

$(BIN): $(F_OBJ)
	$(LDR) $(LDFLAGS) -o $(BIN) $(F_OBJ) -llapack

.F.o:
	$(FC) $(FFLAGS) -c $*.F -llapack

clean:
	rm -f $(F_OBJ) $(C_OBJ) $(BIN) core *.i *.lst
