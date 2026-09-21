#!/usr/bin/env python3
"""Render position CSVs with VMD/Tachyon and FFmpeg; Python 3.8+, no pip packages.

Usage:
  module load vmd/1.9.3 ffmpeg/5.0.1
  python3 PlotDirectory2.py DATA --limit 100 --output preview.mp4
  python3 PlotDirectory2.py DATA
  python3 PlotDirectory2.py DATA --ori-template 'ori_{frame}.csv'

Position files default to pos_*.csv, sorted by their final integer. Each file
has no header, exactly three columns: first row sets the slab z, remaining
rows are particle coordinates. The slab spans [0,20] x [0,20], thickness 1;
particle radius is 0.5. The fixed isometric view fits the FIRST selected frame.

Without --ori-template all particles are red. With it, each matching ori file
must contain exactly one numeric label per particle (flattened CSV). Template
fields: {frame} = final integer preserving leading zeros, {stem}, {name}.
Use --ori-template ori.csv for a constant label file. Number of types comes
from g.csv[0,0], or --types. Colour thresholds reproduce Mathematica: <5000
red, <10000 blue, <15000 purple, otherwise orange, capped at the last type.

--stride subsamples the input; --limit caps the number of output frames.
--threads must not exceed your scheduler's CPU allocation (default 1).
VMD persists across frames. Only one temporary TGA and one Tcl scene exist;
RGB frames go directly to FFmpeg. A successful run also saves the first frame
as OUTPUT.first-frame.ppm and a VMD log as OUTPUT.vmd.log. PPM opens in normal
image viewers, or convert with ffmpeg -i OUTPUT.first-frame.ppm preview.png.
Existing movies are protected unless --overwrite is passed. Failed runs keep
the log but remove the incomplete movie. Run each simulation directory as a
separate scheduler job. No X, GPU, NumPy, PyVista, or Conda is required.
"""

import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import uuid


def natural_key(path):
    nums = re.findall(r"\d+", path.name)
    return (int(nums[-1]) if nums else -1, path.name)


def read_positions(path):
    rows = []
    with path.open(newline="") as f:
        for line, row in enumerate(csv.reader(f), 1):
            if not row or all(not s.strip() for s in row):
                continue
            if len(row) != 3:
                raise ValueError("{}:{}: expected 3 columns".format(path, line))
            values = tuple(float(s) for s in row)
            if not all(math.isfinite(v) for v in values):
                raise ValueError("{}:{}: nonfinite coordinate".format(path, line))
            rows.append(values)
    if not rows:
        raise ValueError("Empty position file: {}".format(path))
    return rows[0][2], rows[1:]


def read_labels(path, count):
    with path.open(newline="") as f:
        labels = [float(s) for row in csv.reader(f) for s in row if s.strip()]
    if len(labels) != count or not all(math.isfinite(v) for v in labels):
        raise ValueError("{}: expected {} finite particle labels, got {}".format(
            path, count, len(labels)))
    return labels


def colour(label, types):
    index = 0 if label < 5000 else 1 if label < 10000 else 2 if label < 15000 else 3
    return (1, 0, 11, 3)[min(index, types - 1)]  # VMD red, blue, purple, orange


def tcl_string(value):
    s = str(value)
    for old, new in (("\\", "\\\\"), ('"', '\\"'), ("$", "\\$"),
                     ("[", "\\["), ("]", "\\]"), ("\n", "\\n"), ("\r", "\\r")):
        s = s.replace(old, new)
    return '"' + s + '"'


def view_point(p):
    x, y, z = p
    return ((x-y)/math.sqrt(2), (-x-y+2*z)/math.sqrt(6), (x+y+z)/math.sqrt(3))


def vertex(p):
    return "{" + " ".join("{:.9g}".format(v) for v in view_point(p)) + "}"


def scene_text(slab_z, points, colours, args, image, first):
    lines = ["graphics $particles delete all", "graphics $slab delete all"]
    last = None
    for p, c in zip(points, colours):
        if c != last:
            lines.append("graphics $particles color {}".format(c))
            last = c
        lines.append("graphics $particles sphere {} radius {} resolution 12".format(
            vertex(p), args.radius))
    b, low, high = args.box_size, slab_z-args.slab_thickness/2, slab_z+args.slab_thickness/2
    corners = [(0,0,low), (b,0,low), (b,b,low), (0,b,low),
               (0,0,high), (b,0,high), (b,b,high), (0,b,high)]
    # Outward winding, two triangles per face; slab is a separate material.
    faces = [(0,3,2,1), (4,5,6,7), (0,1,5,4),
             (1,2,6,5), (2,3,7,6), (3,0,4,7)]
    lines.append("graphics $slab color silver")
    for a, b, c, d in faces:
        for tri in ((a,b,c), (a,c,d)):
            lines.append("graphics $slab triangle " + " ".join(vertex(corners[i]) for i in tri))
    if first:
        lines.extend(["display resetview", "scale by 0.85"])
    # Empty post-render command prevents VMD from launching an image viewer.
    lines.append("render TachyonInternal {} {{}}".format(tcl_string(image)))
    return "\n".join(lines) + "\n"


def read_tga(path, width, height):
    """Read uncompressed/RLE true-colour TGA to top-down RGB, without Pillow."""
    data = path.read_bytes()
    if len(data) < 18:
        raise ValueError("Renderer produced a truncated TGA")
    idlen, cmap, kind = data[:3]
    w, h, depth, flags = struct.unpack_from("<HHBB", data, 12)
    if cmap or kind not in (2, 10) or depth not in (24, 32):
        raise ValueError("Unsupported TGA format: type={}, depth={}".format(kind, depth))
    if (w, h) != (width, height):
        raise ValueError("VMD rendered {}x{}, requested {}x{}".format(w,h,width,height))
    step, count, offset = depth//8, w*h, 18+idlen
    size = count*step
    if kind == 2:
        raw = data[offset:offset+size]
    else:
        raw = bytearray()
        while len(raw) < size:
            header = data[offset]
            offset += 1
            n = (header & 127)+1
            if header & 128:
                raw.extend(data[offset:offset+step]*n)
                offset += step
            else:
                raw.extend(data[offset:offset+n*step])
                offset += n*step
    if len(raw) != size:
        raise ValueError("Truncated or invalid TGA pixel data")
    rgb = bytearray(count*3)
    rgb[0::3], rgb[1::3], rgb[2::3] = raw[2::step], raw[1::step], raw[0::step]
    rows = [rgb[i*w*3:(i+1)*w*3] for i in range(h)]
    if not flags & 32:
        rows.reverse()
    if flags & 16:
        rows = [b"".join(row[i:i+3] for i in range(len(row)-3,-1,-3)) for row in rows]
    return b"".join(rows)


class VMD:
    def __init__(self, args, folder, log):
        self.folder, self.log = folder, log
        env = dict(os.environ, VMDFORCECPUCOUNT=str(args.threads))
        self.proc = subprocess.Popen(
            ["vmd", "-dispdev", "text", "-size", str(args.width), str(args.height),
             "-startup", os.devnull], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1, env=env)

    def run(self, code):
        source = self.folder / "scene.tcl"
        source.write_text(code)
        marker = "PYVMD_" + uuid.uuid4().hex
        command = ('if {[catch {source %s} err]} {puts "%s ERROR $err"} '
                   'else {puts "%s OK"}; flush stdout\n') % (tcl_string(source), marker, marker)
        self.proc.stdin.write(command)
        self.proc.stdin.flush()
        while True:
            line = self.proc.stdout.readline()
            if not line:
                raise RuntimeError("VMD exited unexpectedly; see the .vmd.log file")
            self.log.write(line)
            if marker in line:
                self.log.flush()
                result = line.split(marker, 1)[1].strip()
                if result != "OK":
                    raise RuntimeError("VMD: " + result)
                return

    def close(self):
        if self.proc.poll() is None:
            try:
                self.proc.stdin.write("quit\n")
                self.proc.stdin.flush()
                self.proc.communicate(timeout=10)
            except (OSError, subprocess.TimeoutExpired):
                self.proc.kill()
                self.proc.communicate()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--pattern", default="pos_*.csv")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--fps", type=float, default=30)
    parser.add_argument("--width", type=int, default=1280)
    parser.add_argument("--height", type=int, default=720)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--radius", type=float, default=0.5)
    parser.add_argument("--box-size", type=float, default=20)
    parser.add_argument("--slab-thickness", type=float, default=1)
    parser.add_argument("--ori-template")
    parser.add_argument("--types", type=int, choices=(1,2,3,4))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    for name in ("stride", "fps", "width", "height", "threads", "radius", "box_size", "slab_thickness"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            parser.error("{} must be positive and finite".format(name))
    if args.width % 2 or args.height % 2:
        parser.error("width and height must be even for H.264")
    if args.limit is not None and args.limit < 1:
        parser.error("limit must be positive")
    for program in ("vmd", "ffmpeg"):
        if not shutil.which(program):
            parser.error("{} not found; load the cluster modules first".format(program))
    folder = args.folder.resolve()
    files = sorted((p for p in folder.glob(args.pattern) if p.is_file()), key=natural_key)[::args.stride]
    if args.limit:
        files = files[:args.limit]
    if not files:
        parser.error("No files matching {} in {}".format(args.pattern, folder))
    output = (args.output or folder / "output.mp4").resolve()
    if output.exists() and not args.overwrite:
        parser.error("{} already exists; choose another --output or use --overwrite".format(output))
    if not output.parent.is_dir():
        parser.error("Output directory does not exist: {}".format(output.parent))
    types = args.types
    if args.ori_template and types is None:
        with (folder / "g.csv").open(newline="") as f:
            value = float(next(csv.reader(f))[0])
        if value not in (1,2,3,4):
            parser.error("g.csv[0,0] must be 1, 2, 3, or 4")
        types = int(value)
    preview = output.with_suffix(".first-frame.ppm")
    logpath = output.with_suffix(".vmd.log")
    print("Rendering {} frames -> {}".format(len(files), output), flush=True)
    if not args.ori_template:
        print("No --ori-template: all particles will be red.", flush=True)
    start = time.monotonic()
    # Stage the movie on the output filesystem; replace it only after success.
    fd, partname = tempfile.mkstemp(prefix=output.stem+"-", suffix=".part.mp4", dir=output.parent)
    os.close(fd)
    part = Path(partname)
    vmd = encoder = None
    try:
        with tempfile.TemporaryDirectory(prefix="vmd-render-") as tmp, logpath.open("w") as log:
            tmp = Path(tmp)
            vmd = VMD(args, tmp, log)
            vmd.run('''
display projection orthographic
display depthcue off
display shadows off
display ambientocclusion off
axes location off
color Display Background white
color change rgb 1 1 0 0
color change rgb 0 0 0 1
color change rgb 11 0.5 0 0.5
color change rgb 3 1 0.5 0
set particles [mol new]
set slab [mol new]
material add ParticleMaterial
material change opacity ParticleMaterial 1.0
material change ambient ParticleMaterial 0.3
material change diffuse ParticleMaterial 0.7
material add SlabMaterial
material change opacity SlabMaterial 0.5
material change ambient SlabMaterial 0.3
material change diffuse SlabMaterial 0.7
graphics $particles material ParticleMaterial
graphics $slab material SlabMaterial
''')
            encoder = subprocess.Popen([
                "ffmpeg", "-hide_banner", "-loglevel", "warning", "-y",
                "-f", "rawvideo", "-pixel_format", "rgb24", "-video_size",
                "{}x{}".format(args.width,args.height), "-framerate", str(args.fps),
                "-i", "pipe:0", "-an", "-c:v", "libx264", "-preset", "fast",
                "-crf", "20", "-threads", str(args.threads), "-pix_fmt", "yuv420p",
                "-movflags", "+faststart", str(part)], stdin=subprocess.PIPE)
            image = tmp / "frame.tga"
            for i, path in enumerate(files):
                slab_z, points = read_positions(path)
                if args.ori_template:
                    nums = re.findall(r"\d+", path.stem)
                    labelpath = path.parent / args.ori_template.format(
                        frame=nums[-1] if nums else "", stem=path.stem, name=path.name)
                    colours = [colour(v, types) for v in read_labels(labelpath, len(points))]
                else:
                    colours = [1]*len(points)
                if image.exists():
                    image.unlink()
                vmd.run(scene_text(slab_z, points, colours, args, image, i == 0))
                rgb = read_tga(image, args.width, args.height)
                if i == 0:
                    preview.write_bytes("P6\n{} {}\n255\n".format(args.width,args.height).encode()+rgb)
                encoder.stdin.write(rgb)
                if i == 0 or (i+1) % 10 == 0 or i+1 == len(files):
                    elapsed = time.monotonic()-start
                    print("{}/{} frames, {:.2f} frames/s, elapsed {:.1f}s".format(
                        i+1,len(files),(i+1)/elapsed,elapsed), flush=True)
            encoder.stdin.close()
            if encoder.wait() != 0:
                raise RuntimeError("FFmpeg failed; see its output above")
            vmd.close()
            vmd = None
        if output.exists() and not args.overwrite:
            raise FileExistsError("Output appeared during rendering: {}".format(output))
        os.replace(part, output)
        print("Saved {}\nFirst frame: {}\nVMD log: {}".format(output,preview,logpath))
    finally:
        if vmd is not None:
            vmd.close()
        if encoder is not None and encoder.poll() is None:
            encoder.terminate()
            try:
                encoder.wait(timeout=10)
            except subprocess.TimeoutExpired:
                encoder.kill()
                encoder.wait()
        if part.exists():
            part.unlink()


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as exc:
        print("ERROR: {}".format(exc or "Interrupted"), file=sys.stderr)
        sys.exit(1)
