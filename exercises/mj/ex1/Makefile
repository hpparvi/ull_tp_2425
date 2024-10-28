FC = gfortran

FCFLAGS = -Wall -O2 -g

MODULES = geometry.f90 particle.f90
PROGRAM = ex1.f90

OBJECTS = $(MODULES:.f90=.o) $(PROGRAM:.f90=.o)

EXECUTABLE = leapfrog

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

.PHONY: all clean
