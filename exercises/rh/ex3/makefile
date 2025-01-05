
FC = mpifort
FCFLAGS = -c
OBJECTS = geometry.o particle.o types_mpi.o tree.o e3.o


all: e3

e3: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)
    
%.o: %.f90
	$(FC) $(FCFLAGS) $<
    
clean:
	rm -f *.o *.mod e3
    
test: e3
	./e3


    
