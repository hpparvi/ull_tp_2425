# Compiler
FC = mpif90

# Source and Object Files
SRC = src/geometry.f90 src/particle.f90 src/i_o_utils.f90 src/octree.f90 src/mpi_utils.f90 src/ex3.f90
OBJ = obj/geometry.o obj/particle.o obj/i_o_utils.o obj/octree.o obj/mpi_utils.o obj/ex3.o

# Executable
EXEC = bin/ex3

# Compilation and Linking Flags
CFLAGS = -fdefault-real-8  # Double precision
LFLAGS = -lmpi  # linking mpi flags                   

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

obj/mpi_utils.o: src/mpi_utils.f90 obj/particle.o obj/geometry.o
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/i_o_utils.o: src/i_o_utils.f90 obj/particle.o obj/geometry.o obj/mpi_utils.o 
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/octree.o: src/octree.f90 obj/particle.o obj/geometry.o obj/mpi_utils.o 
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@

obj/ex3.o: src/ex3.f90 obj/geometry.o obj/particle.o obj/i_o_utils.o obj/octree.o obj/mpi_utils.o 
	mkdir -p obj
	$(FC) $(CFLAGS) -c $< -o $@


NP ?= 8
CONFIG ?= ../config_files/config_uniform_collapse.txt

run: $(EXEC)
	mpirun -np $(NP) ./$(EXEC) $(CONFIG)

clean:
	rm -rf obj/*.o bin/ex3 *.mod
