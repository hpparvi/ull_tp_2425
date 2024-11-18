FC = gfortran
FCFLAGS = -c -Jbin
SRC_DIR = src
SRC = $(SRC_DIR)/geometry.f90 $(SRC_DIR)/particles.f90 $(SRC_DIR)/init_conds.f90 $(SRC_DIR)/read_par.f90 $(SRC_DIR)/ex2.f90
OBJ_DIR = bin
OBJECTS = $(OBJ_DIR)/geometry.o $(OBJ_DIR)/particles.o $(OBJ_DIR)/init_conds.o $(OBJ_DIR)/read_par.o $(OBJ_DIR)/ex2.o
PROG := ex2

# to have double precision on all real values
#SINGLE_PRECISION_FLAG =
DOUBLE_PRECISION_FLAG = -fdefault-real-8

all: $(PROG)

test: $(PROG)
	./$(PROG)

$(PROG): $(OBJECTS)
	$(FC) $(DOUBLE_PRECISION_FLAG) -o $@ $(OBJECTS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FCFLAGS) $(DOUBLE_PRECISION_FLAG) -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod

allclean:
	rm -f $(PROG)
