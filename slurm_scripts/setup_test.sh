#!/bin/bash

#SBATCH --job-name=setup-test
#SBATCH --output=setup_test.out
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1

#SBATCH --chdir=..

julia --version
julia --project=. 1_setup.jl -t auto
