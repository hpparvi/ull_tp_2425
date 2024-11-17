# Makefile para el programa ex1

# Compilador de Fortran
FC = gfortran
# Opciones de compilación
FCFLAGS = -O2 -Wall

# Archivos fuente
SRC = ex1.f90 geometry.f90 particle.f90

# Archivos objeto
OBJ = $(SRC:.f90=.o)

# Nombre del ejecutable
EXEC = ex1

# Regla por defecto
all: $(EXEC)

# Regla para enlazar el ejecutable
$(EXEC): $(OBJ)
	$(FC) $(OBJ) -o $(EXEC)

# Regla para compilar los archivos objeto
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Limpiar archivos objeto y ejecutables
clean:
	rm -f $(OBJ) $(EXEC)

.PHONY: all clean
