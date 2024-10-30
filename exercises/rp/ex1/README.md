# Direct Integration N-Body in Fortran

N-body simulator written in Fortran. Initial conditions are either read from file or provided in terminal. Output is saved. Optionally, plot the orbits with a python script. 

## Structure 

Fortran:
- `bin/`: Contains the executable 
- `obj/`: Contains the object files 
- `src/`: Program and modules
  - `src/ex1.f90`: Main program
  - `src/geometry.f90`: Module for vector operations
  - `src/particle.f90`: Module for type particle
  - `src/ics_module.f90`: Module for reading/generating ics

Initial Conditions and Postprocess:
- `python/`: Python files for generating ICS and plotting orbits.
- `data/ics`: Files with initial conditions
- `data/results`: Output data (generated during simulation)
- `figures/broucke_a2`: Example figure of Broucke A2 orbit

## Example Usage 
We use Broucke A2 three-body periodic solution as example for our simulation. The ICS are retrieved from this [webpage](https://observablehq.com/@rreusser/periodic-planar-three-body-orbits)

1. Generate the ICS running `python/make_ics_broucke_a2.py`
2. Compile with `make`
3. Execute `./bin/ex1` and follow prompt instructions
4. Plot and save data with `python/plot_broucke_a2.py`

## Some conlusions
We know the integrator is working because the Broucke A2 three-body orbit looks as it should. 
