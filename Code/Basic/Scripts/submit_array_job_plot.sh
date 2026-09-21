#!/bin/bash
#$ -cwd
#$ -N render_movies
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
#$ -l h_rt=01:59:59,h_data=2G
#$ -t 1-30:1

# Save this file and PlotDirectory2.py in ~/Chemistry/Code/Basic/Scripts/.
# Submit from any directory; SGE writes job logs in the submission directory.
# Full array: qsub ~/Chemistry/Code/Basic/Scripts/submit_job.sh
# Preview: qsub -t 1 -v LIMIT=100 ~/Chemistry/Code/Basic/Scripts/submit_job.sh
# Every frame: qsub -v STRIDE=1 ~/Chemistry/Code/Basic/Scripts/submit_job.sh
# Colour labels default to div_i={frame}.csv; number of types comes from g.csv.
# ORI_TEMPLATE can override the label filename template if needed.
# Default stride 10 covers the whole trajectory. Preview jobs use a separate
# output name, so a preview will not cause the full movie to be skipped.
# One CPU per task: no parallel environment is requested, renderer uses 1 thread.

set -eo pipefail

on_exit() {
    status=$?
    echo "Job ${JOB_ID:-local}, task ${SGE_TASK_ID:-unset} ended: $(date), status=$status"
}
trap on_exit EXIT
echo "Job ${JOB_ID:-local}, task ${SGE_TASK_ID:-unset} started on $(hostname -s): $(date)"

. /u/local/Modules/default/init/modules.sh
module load gcc/11.3.0
module load python/3.9.6
module load vmd/1.9.3
module load ffmpeg/5.0.1
set -u

base_dir="/u/scratch/d/dinoo/GeneticTry19"
script="$HOME/Chemistry/Code/Basic/Scripts/PlotDirectory2.py"
data_dir="${base_dir}/den${SGE_TASK_ID:?Submit as an SGE array job}"
stride="${STRIDE:-10}"
limit="${LIMIT:-}"
ori_template='div_i={frame}.csv'
ori_template="${ORI_TEMPLATE-$ori_template}"
width="${WIDTH:-1280}"
height="${HEIGHT:-720}"

[[ -f "$script" ]] || { echo "Missing script: $script" >&2; exit 1; }
[[ -d "$data_dir" ]] || { echo "Missing directory: $data_dir" >&2; exit 1; }
cd "$data_dir"

output="movie_4views.mp4"
if [[ -n "$limit" ]]; then
    [[ "$limit" =~ ^[1-9][0-9]*$ ]] || { echo "LIMIT must be a positive integer" >&2; exit 1; }
    output="movie_4views_preview_${limit}.mp4"
fi
if [[ -s "$output" ]]; then
    echo "Already exists, skipping: $data_dir/$output"
    exit 0
fi

# Python and its CSV/struct modules require no user-installed packages.
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Use the scheduler's temporary directory for the repeatedly reused TGA/Tcl.
# If it did not supply one, Python chooses its usual temporary location.
if [[ -n "${TMPDIR:-}" && ! -d "$TMPDIR" ]]; then
    mkdir -p "$TMPDIR"
fi

args=("$script" "$data_dir" --pattern 'pos_i=*.csv' --layout quad
      --stride "$stride" --fps 30 --width "$width" --height "$height"
      --threads 1 --output "$output")
if [[ -n "$limit" ]]; then args+=(--limit "$limit"); fi
if [[ -n "$ori_template" ]]; then args+=(--ori-template "$ori_template"); fi
if [[ -n "${TYPES:-}" ]]; then args+=(--types "$TYPES"); fi

echo "Rendering $data_dir -> $output (stride $stride, ${width}x${height})"
/usr/bin/time -v python3 "${args[@]}"
