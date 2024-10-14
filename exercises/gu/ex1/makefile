FC = gfortran
FCFLAGS = -c
OBJECTS : geometry.o particles.o test.o

all: test

test : test

test : geometry.o particles.o test.o
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm *.o *.mod

test : test
	./test
