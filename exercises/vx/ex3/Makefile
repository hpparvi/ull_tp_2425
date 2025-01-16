FC = mpifort
FCFLAGS = -c
OBJECTS = geometry.o definitions.o datatypes_mpi.o env_mpi.o barnes_hut.o main.o
NP = 2

all: main

main: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

run_main : main
	mpirun -np $(NP) ./main

clean:
	rm -f *.o *.mod ex3