# Compiler
FC = gfortran

# Source and Object Files
SRC = src/geometry.f90 src/particle.f90 src/i_o_utils.f90 src/octree.f90 src/ex2.f90
OBJ = obj/geometry.o obj/particle.o obj/i_o_utils.o obj/octree.o obj/ex2.o

# Executable
EXEC = bin/ex2

# Compilation and Linking Flags
CFLAGS = -fopenmp -fdefault-real-8  # Double precision + OpenMP flags
LFLAGS = -fopenmp                   # OpenMP for linking

# Targets
all: $(EXEC)

$(EXEC): $(OBJ)
	mkdir -p bin
	$(FC) $(LFLAGS) -o $(EXEC) $(OBJ)

obj/geometry.o: src/geometry.f90
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/particle.o: src/particle.f90 obj/geometry.o
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/i_o_utils.o: src/i_o_utils.f90 obj/particle.o obj/geometry.o
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/octree.o: src/octree.f90 obj/particle.o obj/geometry.o
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/ex2.o: src/ex2.f90 obj/geometry.o obj/particle.o obj/i_o_utils.o obj/octree.o
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf obj/*.o bin/ex2 *.mod
