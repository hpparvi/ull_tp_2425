from PIL import Image
import subprocess
import os
import sys


class ImageToVideoConverter:
    @staticmethod
    def png_to_mp4(
        fold,
        title="video",
        fps=36,
        digit_format="04d",
        res=None,
        resize_factor=1,
        custom_bitrate=None,
        extension=".jpg",
        dic_guardado='./',
    ):
        files = [f for f in os.listdir(fold) if f.endswith(extension)]
        files.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]))
        print(files)

        if not files:
            raise ValueError("No PNG files found in the specified folder.")

        im = Image.open(os.path.join(fold, files[0]))
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
        output_file = dic_guardado + title + ".mp4"

        crf = 1  # Lower for higher quality, higher for lower quality
        bitrate = custom_bitrate if custom_bitrate else "5000k"
        preset = "slow"
        tune = "film"

        command = f'{ffmpeg_path} -y -r {fps} -i {os.path.join(fold, f"{basename}_%{digit_format}{extension}")} -c:v libx264 -profile:v high -crf {crf} -preset {preset} -tune {tune} -b:v {bitrate} -pix_fmt yuv420p -vf scale={resx}:{resy} {output_file}'

        try:
            result = subprocess.run(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                check=True,
            )
            print('Done!')
        except subprocess.CalledProcessError as e:
            print("Error during video conversion:", e)


if __name__ == "__main__":

    if len(sys.argv) != 6:
        print(
            "Use: python image_to_video.py directory_images/ fps directory_video title resize"
        )
        sys.exit(1)

    dic = str(sys.argv[1])
    fps = sys.argv[2]
    dic_guardado = str(sys.argv[3])
    titulo = sys.argv[4]
    resize = float(sys.argv[5])

    if not (os.path.exists(dic_guardado)):
        os.mkdir(dic_guardado)

    ImageToVideoConverter.png_to_mp4(
        dic,
        extension=".png",
        digit_format="01d",
        fps=fps,
        title=titulo,
        resize_factor=resize,
        dic_guardado=dic_guardado,
    )
