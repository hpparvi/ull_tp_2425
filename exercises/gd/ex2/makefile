FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o ex2_deriv_types.o barnes_hut.o ex2.o
PARTICULAS = geometry.o ex2_deriv_types.o particulas.o

all: ex2 particulas

particulas: $(PARTICULAS)
	$(FC) -o $@ $(PARTICULAS)

ex2 : $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm *.o *.mod ex2 particulas
clean_all :
	rm *.o *.mod *.dat ex2 particulas

test : ex2
	./ex2
	
init: particulas
	./particulas
