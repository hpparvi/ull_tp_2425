FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o input_data.o calcs.o barnes.o ex2.o

all: input_data ex2

ex2: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

input_data: geometry.o particle.o input_data.o
	$(FC) -o $@ geometry.o particle.o input_data.o
%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean: 
	del /Q *.o *.mod

test: input_data ex2
	input_data
	ex2
