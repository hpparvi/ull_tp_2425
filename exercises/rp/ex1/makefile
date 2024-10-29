FC = gfortran  
SRC = src/geometry.f90 src/particle.f90 src/ics_module.f90 src/ex1.f90  
OBJ = obj/geometry.o obj/particle.o obj/ics_module.o obj/ex1.o  
EXEC = bin/ex1  

all: $(EXEC)

$(EXEC): $(OBJ)
	mkdir -p bin
	$(FC) -o $(EXEC) $(OBJ)


obj/geometry.o: src/geometry.f90
	mkdir -p obj
	$(FC) -c src/geometry.f90 -o obj/geometry.o


obj/particle.o: src/particle.f90 obj/geometry.o
	mkdir -p obj
	$(FC) -c src/particle.f90 -o obj/particle.o


obj/ics_module.o: src/ics_module.f90
	mkdir -p obj
	$(FC) -c src/ics_module.f90 -o obj/ics_module.o


obj/ex1.o: src/ex1.f90 obj/geometry.o obj/particle.o obj/ics_module.o
	mkdir -p obj
	$(FC) -c src/ex1.f90 -o obj/ex1.o


clean:
	rm -rf obj/*.o bin/ex1 geometry.mod particle.mod ics_module.mod
