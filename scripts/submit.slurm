#!/bin/bash

#SBATCH --time=12:00:00             # Time limit
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Tasks per node
#SBATCH --cpus-per-task=64          # CPUs per task (adjust as needed)
#SBATCH --mem=256G                  # Memory

julia --project=~/ltmp_calibration --threads=64 "$@"
