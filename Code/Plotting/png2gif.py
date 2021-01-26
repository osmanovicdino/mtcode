from PIL import Image
import glob

# Create the frames
frames = []
imgs = glob.glob("frames/*.png")
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)

# Save into a GIF file that loops forever
frames[0].save('out1.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=125, loop=0)
