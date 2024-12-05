FC = gfortran
FCFLAGS = -c

OBJECTS = geometry.o definitions.o barnes_hut.o main.o


all: main

main: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	del /Q *.o *.mod #it would be rm -f on linux but on windows it's del /Q

test: main
	main