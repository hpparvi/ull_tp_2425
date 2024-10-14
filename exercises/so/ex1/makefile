FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particles.o
PROG := ex1

all: $(PROG)

test: $(PROG)
	./$(PROG)

$(PROG): $(OBJECTS)

#particles: $(OBJECTS)
#	$(FC) -o particles $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm -f *.o *.mod

allclean:
	rm -f particles *.o *.mod
