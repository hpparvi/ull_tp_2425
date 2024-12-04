FC = gfortran
FCFLAGS = -c
LDFLAGS = 
OBJECTS = geometry.o ex2_deriv_types.o barnes_hut.o ex2.o
PARTICULAS = geometry.o ex2_deriv_types.o particulas.o

all: ex2 particulas

particulas: $(PARTICULAS)
	$(FC) -o $@ $(PARTICULAS)

ex2 : $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

ex2_pall: FCFLAGS += -fopenmp
ex2_pall: LDFLAGS += -fopenmp
ex2_pall: $(OBJECTS)
	$(FC) $(LDFLAGS) -o ex2 $(OBJECTS)
	
%.o : %.f90
	$(FC) $(FCFLAGS) $<
	
	
clean :
	rm *.o *.mod ex2 particulas
clean_all :
	rm *.o *.mod *.dat ex2 particulas

test : ex2
	./ex2
	
test_pall: ex2_pall
		
init: particulas
	./particulas
