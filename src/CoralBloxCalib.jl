module CoralBloxCalib

using ADRIA
using DataFrames

include("common/constants.jl")
include("common/types.jl")

include("plot/plotting.jl")

export LocationDataStore,
       ltmp_cover_idx_to_domain,
       composition_idx_to_domain,
       get_ltmp_loc_unique_id,
       get_composition_loc_unique_id,
       extract_param_group_idx,
       extract_param_group

end
