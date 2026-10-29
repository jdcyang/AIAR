FC = ifx
FCFLAGS = -O0 -g -implicitnone -qmkl=sequential -traceback -check all -fpe0 -warn all
LDFLAGS = -O0 -g -qmkl=sequential

OUT = aiar_dbg

MOD = Parameters.o Globals.o

SUBR = SetParameters.o Grids.o TaskSpace.o IterateBellman.o HJBUpdate.o \
       StationaryDistribution.o ComputeGinis.o SolveSteadyState.o SaveOutput.o

OBJ = $(MOD) $(SUBR)

$(OUT).out: $(OBJ) Main.o
	$(FC) $(LDFLAGS) -o $@ $^

Main.o: $(OBJ)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
