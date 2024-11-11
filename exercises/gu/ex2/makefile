FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o tree.o

all: test1

test1 : $(OBJECTS) test1.o
	$(FC) -o $@ $(OBJECTS) test1.o

particulas : $(OBJECTS) particulas.o
	$(FC) -o $@ $(OBJECTS) particulas.o
	./particulas

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod

test : test1
	./test1
