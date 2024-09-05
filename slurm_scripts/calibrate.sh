#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpuq
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=256GB
#SBATCH --job-name="calib"
#SBATCH --chdir=..
#SBATCH --output=slurm_scripts/calibration.out

julia --version
julia --project=. -p auto --threads auto run.jl
