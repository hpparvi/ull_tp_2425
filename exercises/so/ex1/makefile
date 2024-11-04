FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particles.o ex1.o
PROG := ex1

# to have double precision on all real values
#SINGLE_PRECISION_FLAG =
DOUBLE_PRECISION_FLAG = -fdefault-real-8

all: $(PROG)

test: $(PROG)
	./$(PROG)

$(PROG): $(OBJECTS)
	$(FC) $(DOUBLE_PRECISION_FLAG) -o $@ $(OBJECTS)

$(OBJECTS): %.o: %.f90
	$(FC) $(FCFLAGS) $(DOUBLE_PRECISION_FLAG) -o $@ $<

clean:
	rm -f *.o *.mod

allclean:
	rm -f particles *.o *.mod
