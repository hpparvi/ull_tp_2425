
FC = gfortran -fopenmp
FCFLAGS = -c
OBJECTS = geometry.o particle.o tree.o e2.o


all: e2

e2: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)
    
%.o: %.f90
	$(FC) $(FCFLAGS) $<
    
clean:
	rm -f *.o *.mod e2
    
test: e2
	./e2


    
