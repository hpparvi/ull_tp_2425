import os
import shutil
import time
import subprocess
import PIL

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py


def plot_projections(pos, idx=0, lim=(-1, 1), alpha=0.5):

    fs, ratio = 5, 3
    fig, axs = plt.subplots(1, 3, figsize=(ratio*fs, fs))

    dims = [(0,1),(0,2),(1,2)]

    N = pos.shape[1]
    s = 5e4/N

    for i,ax in enumerate(axs):

        ax.scatter(pos[idx,:, dims[i][0]],pos[idx,:, dims[i][1]], lw=0.0, s=s, c='k', alpha=alpha)

    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[0], lim[1])

    plt.show()

    return 



def save_projections(pos, savefold, lim=(-1, 1), alpha=0.5):

    if os.path.exists(savefold):
        shutil.rmtree(savefold)
    
    os.makedirs(savefold)

    N = pos.shape[1]

    fs, ratio = 7, 3
    fig, axs = plt.subplots(1, 3, figsize=(ratio*fs, fs))
    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[0], lim[1])

    dims = [(0,1),(0,2),(1,2)]
    idxs = list(range(pos.shape[0]))

    N = pos.shape[1]
    s = 5e4/N

    for idx in idxs:
        scs = []
        for i,ax in enumerate(axs):
            sc = ax.scatter(pos[idx,:, dims[i][0]],pos[idx,:, dims[i][1]], lw=0.0, s=s, c='k', alpha=alpha)
            scs.append(sc)

        fig.savefig(savefold+f'plot_{idx:04d}.jpg', bbox_inches='tight')
    
        for sc in scs:
            sc.remove()

    plt.close()

    return 



def images_to_mp4(
    fold,
    title="video",
    fps=36,
    digit_format="04d",
    res=None,
    resize_factor=1,
    custom_bitrate=None,
    extension=".jpg",
    reverse=False,  # Adding reverse parameter with default value False
):

    # Get a list of all image files in the directory with the specified extension
    files = [f for f in os.listdir(fold) if f.endswith(extension)]
    files.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]))

    if not files:
        raise ValueError("No image files found in the specified folder.")

    im = PIL.Image.open(os.path.join(fold, files[0]))
    resx, resy = im.size

    if res is not None:
        resx, resy = res
    else:
        resx = int(resize_factor * resx)
        resy = int(resize_factor * resy)
        resx += resx % 2  # Ensuring even dimensions
        resy += resy % 2

    basename = os.path.splitext(files[0])[0].split("_")[0]

    ffmpeg_path = "ffmpeg"
    abs_path = os.path.abspath(fold)
    parent_folder = os.path.dirname(abs_path) + os.sep
    output_file = os.path.join(parent_folder, f"{title}.mp4")

    crf = 5  # Lower CRF for higher quality, higher for lower quality
    bitrate = custom_bitrate if custom_bitrate else "5000k"
    preset = "slow"
    tune = "film"

    # Construct the ffmpeg command
    command = f'{ffmpeg_path} -y -r {fps} -i {os.path.join(fold, f"{basename}_%{digit_format}{extension}")} -c:v libx264 -profile:v high -crf {crf} -preset {preset} -tune {tune} -b:v {bitrate} -pix_fmt yuv420p -vf "scale={resx}:{resy}'
    if reverse:
        command += ",reverse"  # Appends the reverse filter if reverse is True
    command += f'" {output_file}'

    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error during video conversion:", e)






# def images_to_gif(
#     fold,
#     title="video",
#     outfold=None,
#     fps=24,
#     digit_format="04d",
#     quality=500,
#     max_colors=256,
#     extension=".jpg",
#     reverse=False,
# ):
#     files = [f for f in os.listdir(fold) if f.endswith(extension)]
#     files.sort()

#     if reverse:
#         files = files[::-1]

#     name = os.path.splitext(files[0])[0]
#     basename = name.split("_")[0]

#     ffmpeg_path = "ffmpeg"

#     if outfold is None:
#         abs_path = os.path.abspath(fold)
#         parent_folder = os.path.dirname(abs_path) + "/"
#     else:
#         parent_folder = outfold
#         if not os.path.exists(parent_folder):
#             os.makedirs(parent_folder)

#     output_file = os.path.join(parent_folder, f"{title}.gif")

#     # Create a palette with limited colors
#     palette_file = os.path.join(parent_folder, "palette.png")
#     input_pattern = os.path.join(fold, f"{basename}_%{digit_format}{extension}")

#     # Generate palette
#     palette_command = (
#         f'{ffmpeg_path} -framerate {fps} -i {input_pattern} '
#         f'-vf "fps={fps},scale={quality}:-1:flags=lanczos,palettegen=max_colors={max_colors}" '
#         f'-y {palette_file}'
#     )
#     subprocess.run(palette_command, shell=True)

#     if not os.path.exists(palette_file):
#         raise RuntimeError("Palette generation failed. Check ffmpeg installation and inputs.")

#     # Generate GIF using the palette
#     gif_command = (
#         f'{ffmpeg_path} -framerate {fps} -i {input_pattern} -i {palette_file} '
#         f'-lavfi "fps={fps},scale={quality}:-1:flags=lanczos [x]; [x][1:v] paletteuse=dither=bayer:bayer_scale=5" '
#         f'-y {output_file}'
#     )
#     subprocess.run(gif_command, shell=True)

#     # Delete the temporary palette file
#     os.remove(palette_file)