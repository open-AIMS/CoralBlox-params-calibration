# AIMS HPC job submission guide

Submit jobs like:

```bash
cd scripts

sbatch --job-name=loc_calib_test --chdir ~/ltmp_calibration/src submit.slurm ../scripts/run_loc_calib.jl > ../scripts/calib_output.log
```

A log of outputs from the optimizer will saved in `/scripts/calib_output.log`.
