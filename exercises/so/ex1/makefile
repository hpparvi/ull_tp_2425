FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particles.o ex1.o
PROG := ex1

# to have double precision on all real values
#SINGLE_PRECISION_FLAG =
#DOUBLE_PRECISION_FLAG = -fdefault-real-8

all: $(PROG)

test: $(PROG)
	./$(PROG)

$(PROG): $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

#particles: $(OBJECTS)
#	$(FC) -o particles $(OBJECTS)

$(OBJECTS): %.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod

allclean:
	rm -f particles *.o *.mod
