# AIMS HPC job submission guide

Submit jobs like:

```bash
cd scripts

sbatch --job-name=loc_calib_test --chdir ~/CoralBlox-params-calibration/src submit.slurm ../scripts/run_loc_calib.jl --output ../scripts/calib_output.log
```

A log of outputs from the optimizer will saved in `/scripts/calib_output.log`.
