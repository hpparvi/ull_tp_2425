
# Compile the modules first
gfortran -c particle.f90 geometry.f90

# Compile the main program with the modules
gfortran -o ex1 ex1.f90 particle.o geometry.o

# Run the compiled program
./ex1

# Running the program and outputs

When executing the program, you will be prompted to enter values directly in the terminal. Each value should be entered one at a time, pressing Enter after each input. For example, if you select the option to simulate with two particles, you will need to enter each particle’s mass, coordinates, and initial velocities in the following order:
For Particle 1:
  - Mass
  - X coordinate
  - Y coordinate
  - Z coordinate
  - Velocity in X direction
  - Velocity in Y direction
  - Velocity in Z direction
For Particle 2:
  - Mass
  - X coordinate
  - Y coordinate
  - Z coordinate
  - Velocity in X direction
  - Velocity in Y direction
  - Velocity in Z direction

Each prompt will ask for a single value in the sequence shown above. Once you input the value, press Enter, and then proceed to the next.

The program outputs a file called orbits.txt with the following structure:
  Column 1: Time instant
  Column 2: Particle ID (e.g., 1 or 2 in the case of 2-body problem)
  Column 3: X coordinate
  Column 4: Y coordinate
  Column 5: Z coordinate
