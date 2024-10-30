FC = gfortran
PYTHON = python
FCFLAGS = -c
OBJECTS = ex1.o geometry.o particle.o
PLOT_IMAGES_SCRIPT = visual/plot_images.py
IMAGE_TO_VIDEO_SCRIPT = visual/image_to_video.py
DELETE_IMAGES_SCRIPT = visual/delete_images.py

# Variables that can be set from the command line
SIM_NAME ?= test_5part_stable_orbit
OUTPUT_FILE ?= output/$(SIM_NAME).txt 
INPUT_FILE ?= ics/ic_$(SIM_NAME).txt
NUM_CORES ?= None
DIRECTORY_IMAGES = output/images_$(SIM_NAME)
DIRECTORY_VIDEO ?= videos/
TITLE ?= $(SIM_NAME)
FPS ?= 30
RESIZE ?= 1
DELETE ?= TRUE

all: ex1

ex1: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

ex1.o: ex1.f90
	$(FC) $(FCFLAGS) ex1.f90

geometry.o: geometry.f90
	$(FC) $(FCFLAGS) geometry.f90

particle.o: particle.f90
	$(FC) $(FCFLAGS) particle.f90

# To make a video. It is needed ffpmeg installed in the system
video: $(PLOT_IMAGES_SCRIPT) $(IMAGE_TO_VIDEO_SCRIPT)
	@echo "Generating images..."
	$(PYTHON) $(PLOT_IMAGES_SCRIPT) $(OUTPUT_FILE) $(INPUT_FILE) $(NUM_CORES)
	@echo "Generating video..."
	$(PYTHON) $(IMAGE_TO_VIDEO_SCRIPT) $(DIRECTORY_IMAGES) $(FPS) $(DIRECTORY_VIDEO) $(TITLE) $(RESIZE)
	$(PYTHON) $(DELETE_IMAGES_SCRIPT) $(DIRECTORY_IMAGES) $(DELETE)
	
clean:
	rm -f $(OBJECTS) ex1