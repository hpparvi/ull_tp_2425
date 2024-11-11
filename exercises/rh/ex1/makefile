
FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o e1.o

# do make 
# then make test (in order to test the executable)
# and make clean (remove module files and the executable)

all: e1

e1: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)
    
%.o: %.f90
	$(FC) $(FCFLAGS) $<
    
clean:
	rm -f *.o *.mod e1
    
test: e1
	./e1


    
