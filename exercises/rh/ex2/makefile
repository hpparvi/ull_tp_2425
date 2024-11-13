
FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o tree.o e2.o

# do make 
# then make test (in order to test the executable)
# and make clean (remove module files and the executable)

all: e2

e2: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)
    
%.o: %.f90
	$(FC) $(FCFLAGS) $<
    
clean:
	rm -f *.o *.mod e2
    
test: e2
	./e2


    
