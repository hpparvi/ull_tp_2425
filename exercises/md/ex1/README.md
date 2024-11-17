# Exercise 1: Basic Dinamical Simulation
Done by David Mirabal Betancort

## Usage:

The main code is located in `ex1.f90`, which contains the logic of the simulation. It uses `particle.f90`, a module that includes the `particle` type and some essential subroutines for the simulation. The `geometry.f90` module defines basic operations between vectors and points.

To compile the code, you can use the following command in the main directory (`ex1/`): 
```sh
make
```

You can use the code with input files or by introducing initial conditions in the terminal.

- With input files:
```sh
ex1 ics/ic_template.txt
```

- With terminal input:
```sh
ex1
```

This will generate output files with the position of the particles at each simulation time in `output/template.txt`.

To create a video, you can do it in two ways:
> [!WARNING]
> To generate a video, you must install [*FFmpeg*](https://www.ffmpeg.org/).

### With make:
```sh
make video SIM_NAME=template
``` 
`template` is the name of the simulation, so the input file should be `ics/ic_template.txt`, the output file `output/template.txt`, and the video `videos/template.mp4`.

You can manually change several video parameters such as `INPUT_FILE`, `OUTPUT_FILE`, `NUM_CORES` (defaults to the maximum available), `DIRECTORY_IMAGES` (the directory where images are saved), `DIRECTORY_VIDEO` (the directory where the video is saved), `TITLE`, `FPS`, `RESIZE` (resize factor to reduce the video size), and `DELETE` (whether to delete images created during video creation; defaults to True).
```sh
make video SIM_NAME=template INPUT_FILE=ics/ic_template.txt OUTPUT_FILE=output/template.txt NUM_CORES=8 DIRECTORY_IMAGES=output/images_template/ DIRECTORY_VIDEO=videos/ TITLE=template FPS=30 RESIZE=1
``` 

### Manually: 
You can run each Python script by hand to make the video:
```sh
python visual/plot_images.py output/template.txt ics/ic_template.txt 8
python visual/image_to_video.py output/images_template/ 30 videos/ template 1
```

8 is the number of cores (threads) to use, 30 is fps and 1 is the resize factor.

