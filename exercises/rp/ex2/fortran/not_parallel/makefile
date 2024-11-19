FC = gfortran  
SRC = src/geometry.f90 src/particle.f90 src/i_o_utils.f90 src/octree.f90 src/ex2.f90  
OBJ = obj/geometry.o obj/particle.o obj/i_o_utils.o obj/octree.o obj/ex2.o  
EXEC = bin/ex2  

# SINGLE_PRECISION_FLAG =
DOUBLE_PRECISION_FLAG = -fdefault-real-8

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

obj/i_o_utils.o: src/i_o_utils.f90 obj/particle.o obj/geometry.o
	mkdir -p obj
	$(FC) -c src/i_o_utils.f90 -o obj/i_o_utils.o

obj/octree.o: src/octree.f90 obj/particle.o obj/geometry.o
	mkdir -p obj
	$(FC) -c src/octree.f90 -o obj/octree.o


obj/ex2.o: src/ex2.f90 obj/geometry.o obj/particle.o obj/i_o_utils.o obj/octree.o
	mkdir -p obj
	$(FC) -c src/ex2.f90 -o obj/ex2.o


clean:
	rm -rf obj/*.o bin/ex2 geometry.mod particle.mod i_o_utils.mod octree.mod
