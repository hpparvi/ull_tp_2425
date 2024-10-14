FC = gfortran
FCFLAGS = -c

all: e1

test: e1
	./e1

e1: m1.o e1.o 
	$(FC) -o e1 e1.o m1.o

e1.o: e1.f90
	$(FC) $(FCFLAGS) -c e1.f90

m1.o: m1.f90
	$(FC) $(FCFLAGS) -c m1.f90

clean:
	rm *.o *.mod e1
