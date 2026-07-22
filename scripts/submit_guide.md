# AIMS HPC job submission guide

See setup instructions in the README.md.

If the path to your calibration repo is `~/CoralBlox-params-calibration/` you can submit jobs using:

```bash
cd scripts

sbatch --job-name=loc_calib_test --chdir ~/CoralBlox-params-calibration/ submit.slurm ../scripts/run_loc_calib.jl --output ../scripts/calib_output.log
```

A log of outputs from the optimizer will saved in `/scripts/calib_output.log`.
