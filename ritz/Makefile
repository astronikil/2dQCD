SHELL = /bin/sh

BIN1 = u1eig
BIN2 = tune

FC = mpif77
LDR = mpif77

FFLAGS  = -funroll-loops -O3
LDFLAGS = -funroll-loops -O3

.SUFFIXES:  .f .c .o

F_OBJ1 =  main.f random.f mathtools.f sort.f spectrum.f startritz.f ritz.f diagritz.f GS.f ferm.F addrc.F io.F wloops.F  miscroutines.f ipr.f
F_OBJ2 =  tune.f random.f mathtools.f sort.f spectrum.f startritz.f ritz.f diagritz.f GS.f ferm.F addrc.F io.F wloops.F miscroutines.f

compile: $(BIN1) $(BIN2)

$(BIN1): 
	sed -ie "s/nev=.*/nev=20)/g" global.f
	$(LDR) $(LDFLAGS) -o $(BIN1) $(F_OBJ1)
	@echo DONE EIG

$(BIN2): 
	sed -ie "s/nev=.*/nev=1)/g" global.f
	$(LDR) $(LDFLAGS) -o $(BIN2) $(F_OBJ2)
	@echo DONE TUNE

clean:
	rm -f $(BIN1) $(BIN2) core *.i *.lst
