FC = mpifort
BASE_FCFLAGS = -Wall -Wextra -c -Jbin
BASE_LDFLAGS =
# OPENMP = yes
# OPENMP_FLAGS = -fopenmp
# Conditional Flags
# ifeq ($(OPENMP), yes)
#    FCFLAGS = $(BASE_FCFLAGS) $(OPENMP_FLAGS)
#    LDFLAGS = $(OPENMP_FLAGS)
#else
#    FCFLAGS = $(BASE_FCFLAGS)
#    LDFLAGS =
#endif
LDFLAGS =
FCFLAGS = $(BASE_FCFLAGS)


SRC_DIR = src
SRC = $(SRC_DIR)/geometry.f90 $(SRC_DIR)/particles.f90 $(SRC_DIR)/init_conds.f90 $(SRC_DIR)/read_par.f90 $(SRC_DIR)/barnes_hut.f90 $(SRC_DIR)/ex3.f90
OBJ_DIR = bin
OBJECTS = $(OBJ_DIR)/geometry.o $(OBJ_DIR)/particles.o $(OBJ_DIR)/init_conds.o $(OBJ_DIR)/read_par.o $(OBJ_DIR)/barnes_hut.o $(OBJ_DIR)/ex3.o
PROG := ex3

# to have double precision on all real values
DOUBLE_PRECISION_FLAG = -fdefault-real-8

all: $(PROG)

test: $(PROG)
	./$(PROG)

$(PROG): $(OBJECTS)
	$(FC) $(LDFLAGS) $(DOUBLE_PRECISION_FLAG) -o $@ $(OBJECTS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FCFLAGS) $(DOUBLE_PRECISION_FLAG) -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod

allclean:
	rm -f $(PROG)
