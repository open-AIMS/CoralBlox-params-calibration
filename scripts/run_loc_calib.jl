src_dir = joinpath(expanduser("~"), "ltmp_calibration", "src")
cd(src_dir)
@info pwd()


include(joinpath(src_dir, "1_setup.jl"))
include(joinpath(src_dir, "2b_location_calibration.jl"))
