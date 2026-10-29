# source /opt/intel/oneapi/setvars.sh before make
FC = ifx
FCFLAGS = -O3 -implicitnone -qmkl=sequential -traceback
LDFLAGS = -O3 -qmkl=sequential

OUT = aiar

MOD = Parameters.o Globals.o

SUBR = SetParameters.o Grids.o TaskSpace.o IterateBellman.o HJBUpdate.o \
       StationaryDistribution.o ComputeGinis.o SolveSteadyState.o SaveOutput.o

OBJ = $(MOD) $(SUBR)

$(OUT).out: $(OBJ) Main.o
	$(FC) $(LDFLAGS) -o $@ $^

Main.o: $(OBJ)
$(SUBR): $(MOD)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
