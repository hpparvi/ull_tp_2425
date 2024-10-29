FC = gfortran  # Fortran compiler
SRC = src/geometry.f90 src/particle.f90 src/ics_module.f90 src/ex1.f90  # List of source files
OBJ = obj/geometry.o obj/particle.o obj/ics_module.o obj/ex1.o  # Object files
EXEC = bin/ex1  # The final executable

all: $(EXEC)

# Link object files to create executable
$(EXEC): $(OBJ)
	mkdir -p bin
	$(FC) -o $(EXEC) $(OBJ)

# Rule to compile geometry.f90
obj/geometry.o: src/geometry.f90
	mkdir -p obj
	$(FC) -c src/geometry.f90 -o obj/geometry.o

# Rule to compile particle.f90 (dependent on geometry.mod)
obj/particle.o: src/particle.f90 obj/geometry.o
	mkdir -p obj
	$(FC) -c src/particle.f90 -o obj/particle.o

# Rule to compile ics_module.f90
obj/ics_module.o: src/ics_module.f90
	mkdir -p obj
	$(FC) -c src/ics_module.f90 -o obj/ics_module.o

# Rule to compile ex1.f90 (dependent on geometry, particle, and ics_module)
obj/ex1.o: src/ex1.f90 obj/geometry.o obj/particle.o obj/ics_module.o
	mkdir -p obj
	$(FC) -c src/ex1.f90 -o obj/ex1.o

# Clean rule
clean:
	rm -rf obj/*.o bin/ex1 geometry.mod particle.mod ics_module.mod
