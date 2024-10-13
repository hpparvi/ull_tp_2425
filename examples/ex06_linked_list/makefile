FC = gfortran
FCFLAGS = -c
OBJECTS = slist.o test.o

all: test_slist

test: test_slist
	./test_slist

test_slist: $(OBJECTS)
	$(FC) -o test_slist $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm -f *.o *.mod *~ test_slist
