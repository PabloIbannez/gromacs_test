#!/usr/bin/env bash
#SBATCH --job-name=gmx-water-fixed
#SBATCH --array=0-14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

set -euo pipefail

CASE_DIRS=(
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference/water_5_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/candidate/water_5_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference_double/water_5_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference/water_10_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/candidate/water_10_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference_double/water_10_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference/water_50_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/candidate/water_50_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference_double/water_50_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference/water_100_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/candidate/water_100_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference_double/water_100_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference/water_150_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/candidate/water_150_3_3"
    "/home/pablo/projects/gromacs_test/water_fixed_precision_campaign/reference_double/water_150_3_3"
)

CASE_DIR="${CASE_DIRS[$SLURM_ARRAY_TASK_ID]}"

echo "Running case $SLURM_ARRAY_TASK_ID in $CASE_DIR"
make -C "$CASE_DIR" all
