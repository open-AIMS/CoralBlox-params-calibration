#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=128GB
#SBATCH --job-name="calib"
#SBATCH --chdir=..
#SBATCH --output=slurm_scripts/calibration.out

julia --version
julia --project=. --threads auto run.jl
