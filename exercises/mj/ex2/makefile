# Fortran compiler
FC = gfortran

# Compiler flags
FCFLAGS = -c
LDFLAGS = -o

# Module and program files
MODULE_1 = geometry.f90
MODULE_2 = particle.f90
MODULE_3 = bh.f90
PROGRAM = ex2.f90

# Object files
OBJ_1 = geometry.o
OBJ_2 = particle.o
OBJ_3 = bh.o
PROGRAM_OBJ = ex2.o

# Output executable name
EXECUTABLE = sim

all: $(EXECUTABLE)

# Compile the executable, ensuring particle.o, bh.o are ready before main.o
$(EXECUTABLE): $(OBJ_1) $(OBJ_2) $(OBJ_3) $(PROGRAM_OBJ)
	$(FC) $(OBJ_1) $(OBJ_2) $(OBJ_3) $(PROGRAM_OBJ) $(LDFLAGS) $(EXECUTABLE)

# Compilation rules for each module and the main program
$(OBJ_1): $(MODULE_1)
	$(FC) $(FCFLAGS) $(MODULE_1)

$(OBJ_2): $(MODULE_2)
	$(FC) $(FCFLAGS) $(MODULE_2)

$(OBJ_3): $(MODULE_3) $(OBJ_2)
	$(FC) $(FCFLAGS) $(MODULE_3)

$(PROGRAM_OBJ): $(PROGRAM) $(OBJ_2) $(OBJ_3)
	$(FC) $(FCFLAGS) $(PROGRAM)

test: $(EXECUTABLE)
	./$(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) *.mod *.o
