FC = gfortran
PYTHON = python
MAKE = make
FCFLAGS = -c -g -fcheck=bounds
PARFLAGS = -fopenmp -g -fcheck=bounds
OBJECTS = ex2.o barnes_hut.o calculations.o utils.o geometry.o particle.o
OBJECTS_PARALLEL = ex2_par.o barnes_hut_par.o calculations_par.o utils.o geometry.o particle_par.o
TESTS = tests/test_stable_orbit_file.py
PLOT_IMAGES_SCRIPT = visual/plot_images.py
IMAGE_TO_VIDEO_SCRIPT = visual/image_to_video.py
DELETE_IMAGES_SCRIPT = visual/delete_images.py

# Variables that can be set from the command line
SIM_NAME ?= test_5part_stable_orbit
OUTPUT_FILE ?= output/$(SIM_NAME).dat
INPUT_FILE ?= ics/ic_$(SIM_NAME).txt
NUM_CORES ?= None
DIRECTORY_IMAGES = output/images_$(SIM_NAME)
DIRECTORY_VIDEO ?= videos/
TITLE ?= $(SIM_NAME)
FPS ?= 30
RESIZE ?= 1
DELETE ?= TRUE

all: ex2

parallel: $(OBJECTS_PARALLEL)
	$(FC) $(PARFLAGS) -o ex2 $(OBJECTS_PARALLEL)

ex2: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

ex2.o: ex2.f90
	$(FC) $(FCFLAGS) ex2.f90

ex2_par.o: ex2.f90
	$(FC) $(PARFLAGS) $(FCFLAGS) -o $@ ex2.f90

geometry.o: geometry.f90
	$(FC) $(FCFLAGS) geometry.f90

particle.o: particle.f90
	$(FC) $(FCFLAGS) particle.f90

particle_par.o: particle.f90
	$(FC) $(PARFLAGS) $(FCFLAGS) -o $@ particle.f90

calculations.o: calculations.f90
	$(FC) $(FCFLAGS) calculations.f90

calculations_par.o: calculations.f90
	$(FC) $(PARFLAGS) $(FCFLAGS) -o $@ calculations.f90

barnes_hut.o: barnes_hut.f90
	$(FC) $(FCFLAGS) barnes_hut.f90

barnes_hut_par.o: barnes_hut.f90
	$(FC) $(PARFLAGS) $(FCFLAGS) -o $@ barnes_hut.f90

utils.o: utils.f90
	$(FC) $(FCFLAGS) utils.f90

# To make a video. It is needed ffpmeg installed in the system
video: $(PLOT_IMAGES_SCRIPT) $(IMAGE_TO_VIDEO_SCRIPT)
	@echo "Generating images..."
	$(PYTHON) $(PLOT_IMAGES_SCRIPT) $(OUTPUT_FILE) $(INPUT_FILE) $(NUM_CORES)
	@echo "Generating video..."
	$(PYTHON) $(IMAGE_TO_VIDEO_SCRIPT) $(DIRECTORY_IMAGES) $(FPS) $(DIRECTORY_VIDEO) $(TITLE) $(RESIZE)
	$(PYTHON) $(DELETE_IMAGES_SCRIPT) $(DIRECTORY_IMAGES) $(DELETE)

all_videos: 
	for file in output/*.dat; do \
		SIM_NAME_V=$$(basename $$file .dat); \
		$(MAKE) video OUTPUT_FILE=$$file SIM_NAME=$$SIM_NAME_V; \
	done
	@echo "All videos generated"

all_ics: 
	for file in ics/*.txt; do \
		./ex2 $$file; \
	done

test: 
	$(PYTHON) $(TESTS)

clean:
	rm -f $(OBJECTS) $(OBJECTS_PARALLEL) ex2