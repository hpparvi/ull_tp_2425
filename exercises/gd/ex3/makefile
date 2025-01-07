FC = mpifort
FCFLAGS = -c
OBJECTS = geometry.o ex3_deriv_types.o barnes_hut.o ex3.o

all: ex3


ex3 : $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

	
%.o : %.f90
	$(FC) $(FCFLAGS) $<
	
	
clean :
	rm *.o *.mod ex3
clean_all :
	rm *.o *.mod *.dat ex3

test : ex3
	mpirun -n 12 ./ex3
	
