# Variables
FC = gfortran  # Fortran compiler
SRC_DIR = src  # Source directory
OBJ_DIR = obj  # Object files directory
BIN_DIR = bin  # Executable directory
SRC = $(wildcard $(SRC_DIR)/*.f90)  # Automatically find all .f90 files in src/
OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC))  # Create obj files from src files
EXEC = $(BIN_DIR)/exec  # The final executable

# Default rule
all: $(EXEC)

# Rule to create the executable
$(EXEC): $(OBJ)
	mkdir -p $(BIN_DIR)
	$(FC) -o $(EXEC) $(OBJ)

# Rule to compile .f90 files into .o files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR)
	$(FC) -c $< -o $@

# Clean rule
clean:
	rm -rf $(OBJ_DIR)/*.o $(BIN_DIR)/exec *.mod
