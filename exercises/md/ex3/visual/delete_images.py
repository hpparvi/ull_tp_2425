import os
import sys

# To delete images when the video is done

file = sys.argv[1]
delete = sys.argv[2]
if __name__ == '__main__':
    if delete.lower() in ['true', 't', 'y', 'yes']:
        print(f"Deleting images in {file}")
        for f in os.listdir(file):
            os.remove(os.path.join(file, f))
        os.rmdir(file)
        print("Done!")
    else:
        print("Images are not going to be deleted")
        print("Done!")
