FC = mpifort
FCFLAGS = -c

all: e0 e1 e2 e3 e4 e5

e0: e0.o
	$(FC) -o e0 e0.o

e1: e1.o
	$(FC) -o e1 e1.o

e2: e2.o
	$(FC) -o e2 e2.o

e3: e3.o
	$(FC) -o e3 e3.o

e4: e4.o
	$(FC) -o e4 e4.o

e5: e5.o
	$(FC) -o e5 e5.o

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm -f *.o *.mod

allclean:
	rm -f e? *.o *.mod
