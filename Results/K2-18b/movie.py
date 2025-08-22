import os
import re
import imageio

# Path to the folder containing your images
folder_path = "live_plot/"

# Pattern to match the image files and extract chapter and frame numbers
pattern = re.compile(r'TP_profile_step_(\d+)_(\d+)\.png')

# Collect all images that match the pattern
images_info = []
for filename in os.listdir(folder_path):
    match = pattern.match(filename)
    if match:
        chapter = int(match.group(1))
        frame = int(match.group(2))
        images_info.append((chapter, frame, filename))

# Sort images first by chapter, then by frame
images_info.sort(key=lambda x: (x[0], x[1]))

if not images_info:
    raise ValueError("No images matching pattern found in the folder.")

# Create the movie
output_file = "output_movie.mp4"
fps = 9

with imageio.get_writer(output_file, fps=fps) as writer:
    for _, _, filename in images_info:
        image_path = os.path.join(folder_path, filename)
        image = imageio.imread(image_path)
        writer.append_data(image)

print(f"Movie successfully created as {output_file}")