FC = gfortran
FCFLAGS = -c
SRC = src/geometry.f90 src/particle.f90 src/ex2.f90 src/init_conds.f90
OBJECTS = geometry.o particles.o ex2.o
PROG := ex2

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
