FC = gfortran
FCFLAGS = -c
OBJECTS = m1.o e1.o

all: e1

test: e1
	./e1

e1: $(OBJECTS)
	$(FC) -o e1 $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm -f *.o *.mod

allclean:
	rm -f e1 *.o *.mod
