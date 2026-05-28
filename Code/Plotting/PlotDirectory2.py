import numpy as np
import pyvista as pv
import pandas as pd
import subprocess
import glob
import os
import re
import sys


def natural_key(path):
    """Sort files by the trailing integer in the filename (the frame index)."""
    nums = re.findall(r'\d+', os.path.basename(path))
    return int(nums[-1]) if nums else 0


def read_frame(path):
    """One CSV per frame. Row 0 is the slab reference (its z sets the slab
    height); rows 1.. are particle coordinates. Returns (slab_z, coords)."""
    arr = pd.read_csv(path, header=None, dtype=np.float32).values
    return float(arr[0, 2]), arr[1:]


def make_video(folder, pattern='data_*.csv', fps=30,
               box_size=20.0, slab_thickness=1.0, radius=0.5,
               sphere_color=(120, 120, 120), fixed_camera=True):

    files = sorted(glob.glob(os.path.join(folder, pattern)), key=natural_key)
    n_frames = len(files)
    if n_frames == 0:
        raise FileNotFoundError(f'no files matching {pattern!r} in {folder}')
    print(f'found {n_frames} frames')

    sphere = pv.Sphere(radius=radius, theta_resolution=8, phi_resolution=8)

    pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    pl.set_background('white')
    pl.enable_lightkit()
    pl.disable_anti_aliasing()

    # Fixed camera framed on the first frame + the 20x20 slab footprint.
    # Mathematica's Graphics3D view is static, so this is the closest analogue.
    slab_z0, pts0 = read_frame(files[0])
    if fixed_camera:
        lo = np.minimum(pts0.min(axis=0), [0, 0, slab_z0])
        hi = np.maximum(pts0.max(axis=0), [box_size, box_size, slab_z0])
        center = (lo + hi) / 2
        dist = float(np.linalg.norm(hi - lo)) * 1.2
        direction = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
        pl.camera.position = center + direction * dist
        pl.camera.focal_point = center
        pl.camera.up = (0, 0, 1)

    pipe = subprocess.Popen(
        ['ffmpeg', '-y', '-f', 'rawvideo', '-pixel_format', 'rgb24',
         '-video_size', '1920x1080', '-framerate', str(fps),
         '-i', '-', '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
         os.path.join(folder, 'output.mp4')],
        stdin=subprocess.PIPE
    )

    for i, path in enumerate(files):
        slab_z, pts = read_frame(path)

        # Spheres (variable count per frame is fine — rebuilt each time).
        cloud = pv.PolyData(pts)
        glyphs = cloud.glyph(geom=sphere, scale=False, orient=False)
        pl.add_mesh(glyphs, color=sphere_color, name='particles')

        # Semi-transparent slab: 20x20 in xy, thickness 1 in z, centered on slab_z.
        slab = pv.Box(bounds=(0, box_size, 0, box_size,
                              slab_z - slab_thickness / 2,
                              slab_z + slab_thickness / 2))
        pl.add_mesh(slab, color='lightgray', opacity=0.5, name='slab')

        # Optional: render.py-style centroid tracking instead of a fixed view.
        if not fixed_camera:
            center = pts.mean(axis=0)
            dist = max(pts.std(axis=0).max() * 10, 50.0)
            direction = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
            pl.camera.position = center + direction * dist
            pl.camera.focal_point = center
            pl.camera.up = (0, 0, 1)

        pipe.stdin.write(pl.screenshot(return_img=True).tobytes())

        if i % 100 == 0:
            print(f'{i}/{n_frames}')

    pipe.stdin.close()
    pipe.wait()
    pl.close()
    print('done')


if __name__ == '__main__':
    folder = sys.argv[1]
    make_video(folder)