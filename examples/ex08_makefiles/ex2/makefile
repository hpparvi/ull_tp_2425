FC = gfortran
FCFLAGS = -c
OBJECTS = m1.o e1.o

all: e1

test: e1
	./e1

e1: $(OBJECTS)
	$(FC) -o e1 $(OBJECTS)

e1.o: e1.f90
	$(FC) $(FCFLAGS) $<

m1.o: m1.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm *.o *.mod

allclean:
	rm e1 *.o *.mod
