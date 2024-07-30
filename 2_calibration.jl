include("common.jl")

# Get coral parameter list
coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)
not_constant = coral_params.is_constant .== false

"""
Model simulations end at 2022, so north/central/south obs argument ignores the last entry
which is for 2023.
"""
function obj_func(
    init_values;
    dom=dom,
    north_cover=ltmp_north[ltmp_north_period, [:lower, :response, :upper]],
    central_cover=ltmp_central[ltmp_central_period, [:lower, :response, :upper]],
    south_cover=ltmp_south[ltmp_south_period, [:lower, :response, :upper]],
    coral_params=coral_params[not_constant, :],
    start_year=start_year,
    end_year=end_year,
)
    n_cover_start = length(dom.init_coral_cover)
    gen_init_cover = reshape(init_values[1:n_cover_start], size(dom.init_coral_cover)...)

    loc_k_areas = site_k_area(dom)
    loc_areas = site_area(dom)

    # Check constraints
    total_cover_gen = vec(sum(gen_init_cover, dims=1))

    north_obs = north_cover.response
    central_obs = central_cover.response
    south_obs = south_cover.response

    # If start condition is out of bounds
    start_score = zeros(3)
    total_north = (total_cover_gen[NORTH_MASK] .* loc_k_areas[NORTH_MASK]) ./ loc_areas[NORTH_MASK]
    start_score[1] = abs(mean(total_north) - north_obs[1])
    if !(north_cover.lower[1] .<= mean(total_north) .<= north_cover.upper[1])
        start_score[1] *= 100
        start_score[1] += 1e6
    end

    total_central = (total_cover_gen[CENTRAL_MASK] .* loc_k_areas[CENTRAL_MASK]) ./ loc_areas[CENTRAL_MASK]
    start_score[2] = abs(mean(total_central) - central_obs[1])
    if !(central_cover.lower[1] .<= mean(total_central) .<= central_cover.upper[1])
        start_score[2] *= 100
        start_score[2] += 1e6
    end

    total_south = (total_cover_gen[SOUTH_MASK] .* loc_k_areas[SOUTH_MASK]) ./ loc_areas[SOUTH_MASK]
    start_score[3] = abs(mean(total_south) - south_obs[1])
    if !(south_cover.lower[1] .<= mean(total_south) .<= south_cover.upper[1])
        start_score[3] *= 100
        start_score[3] += 1e6
    end

    if sum(start_score) >= 1e6
        return sum(start_score)
    end

    # If initial cover pass constraint rules, assign them and run
    dom2 = deepcopy(dom)
    dom2.init_coral_cover.data .= gen_init_cover

    scen = ADRIA.param_table(dom2)
    coral_param_values = init_values[n_cover_start+1:end]
    scen[1, coral_params.fieldname] = coral_param_values

    res = nothing
    try
        res = ADRIA.run_model(dom2, scen[1, :])
    catch err
        return sum(start_score) + 5e5
    end

    north_mean_cover = zeros(size(res.raw, 1))
    center_mean_cover = zeros(size(res.raw, 1))
    south_mean_cover = zeros(size(res.raw, 1))

    # Hacky manual specification of observation years to use to compare against LTMP data
    # which misses some years (ignoring 2023, the last entry)
    ref_years = start_year:end_year
    comp_years_north = (ref_years .∈ [ltmp_north[ltmp_north.Year .>= start_year, :Year]])
    comp_years_center = (ref_years .∈ [ltmp_central[ltmp_central.Year .>= start_year, :Year]])
    comp_years_south = (ref_years .∈ [ltmp_south[ltmp_south.Year .>= start_year, :Year]])

    for ts in axes(res.raw, 1)
        rac = vec(sum(res.raw[ts, :, :], dims=1) .* loc_k_areas') ./ loc_areas
        north_mean_cover[ts] = mean(rac[NORTH_MASK])
        center_mean_cover[ts] = mean(rac[CENTRAL_MASK])
        south_mean_cover[ts] = mean(rac[SOUTH_MASK])
    end

    north_mean_cover = north_mean_cover[comp_years_north]
    center_mean_cover = center_mean_cover[comp_years_center]
    south_mean_cover = south_mean_cover[comp_years_south]

    # Determine performance
    north_perf_series = MAEE_series(north_mean_cover, north_obs)
    central_perf_series = MAEE_series(center_mean_cover, central_obs)
    south_perf_series = MAEE_series(south_mean_cover, south_obs)

    north_perf = temporal_variability(north_perf_series)
    central_perf = temporal_variability(central_perf_series)
    south_perf = temporal_variability(south_perf_series)

    # Penalize poor performance at lowest trough
    north_lowest = argmin(north_obs)
    central_lowest = argmin(central_obs)
    south_lowest = argmin(south_obs)
    trough_score = [north_perf_series[north_lowest], central_perf_series[central_lowest], south_perf_series[south_lowest]] .* 100.0

    # Penalize poor performance at end of time series
    end_score = [north_perf_series[end], central_perf_series[end], south_perf_series[end]] .* 100.0

    return sum([north_perf, central_perf, south_perf]) + sum(trough_score) + sum(end_score)
end

# Define parameter space to scan over
search_space = repeat([(0.0, 0.6)], prod(*, size(dom.init_coral_cover)))

not_linear_ext = .!(occursin.("linear_extension", string.(coral_params.fieldname)))
not_diam = .!(occursin.("mean_colony", string.(coral_params.fieldname)))
selected_coral_params = not_constant .& not_linear_ext .& not_diam

coral_params[selected_coral_params, :lower_bound] .= coral_params[selected_coral_params, :val] .* 0.05
coral_params[selected_coral_params, :upper_bound] .= coral_params[selected_coral_params, :val] .* 10.0

coral_params[(.!not_linear_ext .| .!not_diam), :lower_bound] .= coral_params[(.!not_linear_ext .| .!not_diam), :val] .* 0.05

# Maximum linear extension is diameter of size class
coral_params[.!not_linear_ext, :upper_bound] .= coral_params[.!not_diam, :val]

# Diameters are already set to upper bound
coral_params[.!not_diam, :upper_bound] .= coral_params[.!not_diam, :val]

coral_bounds = Tuple.(eachrow(coral_params[not_constant, [:lower_bound, :upper_bound]]))
append!(search_space, coral_bounds)

best_score_file = "Outputs/best_initial_state_w_coral_$(string(today())).dat"
if !@isdefined(best_init_state) && isfile(best_score_file)
    best_init_state = deserialize(best_score_file)
else
    best_init_state = Float64.(vcat(dom.init_coral_cover.data[:], coral_params[not_constant, :val]))
end

@assert all(first.(search_space) .<= best_init_state .<= last.(search_space)) "Initial state is out of bounds"

if !@isdefined(best_init_state) || isnothing(best_init_state)
    # Include additional config if using BorgMOEA
    # Method=:borg_moea,
    # FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
    res = bboptimize(
        obj_func;
        SearchRange=search_space,
        MaxSteps=1_000_000,
        NThreads=Threads.nthreads()-4
    );
elseif !isnothing(best_init_state)
    res = bboptimize(
        obj_func,
        best_init_state;  # provide an initial solution
        SearchRange=search_space,
        MaxSteps=1_000_000,
        NThreads=Threads.nthreads()-4
    );
end

best_fitness(res)
best_init_state = best_candidate(res)
n_cover_start = length(dom.init_coral_cover)
best_init_cover = reshape(best_init_state[1:n_cover_start], size(dom.init_coral_cover))
sum(best_init_cover, dims=1)

serialize(best_score_file, best_init_state)

dom.init_coral_cover .= best_init_cover
coral_param_values = best_init_state[n_cover_start+1:end]
scens = ADRIA.param_table(dom)
scens[1, coral_params[not_constant, :].fieldname] = coral_param_values

rs_raw = ADRIA.run_model(dom, scens[1, :])
s_rac = (dropdims(sum(rs_raw.raw, dims=2), dims=2) .* site_k_area(dom)') ./ site_area(dom)'

ref_years = start_year:end_year
north_mean_cover = zeros(size(rs_raw.raw, 1))
center_mean_cover = zeros(size(rs_raw.raw, 1))
south_mean_cover = zeros(size(rs_raw.raw, 1))
comp_years_north = (ref_years .∈ [ltmp_north[ltmp_north.Year .>= start_year, :Year]])
comp_years_center = (ref_years .∈ [ltmp_central[ltmp_central.Year .>= start_year, :Year]])
comp_years_south = (ref_years .∈ [ltmp_south[ltmp_south.Year .>= start_year, :Year]])

for ts in axes(rs_raw.raw, 1)
    rac = vec(sum(rs_raw.raw[ts, :, :], dims=1) .* site_k_area(dom)') ./ site_area(dom)
    north_mean_cover[ts] = mean(rac[NORTH_MASK])
    center_mean_cover[ts] = mean(rac[CENTRAL_MASK])
    south_mean_cover[ts] = mean(rac[SOUTH_MASK])
end

north_mean_cover = north_mean_cover[comp_years_north]
center_mean_cover = center_mean_cover[comp_years_center]
south_mean_cover = south_mean_cover[comp_years_south]

north_res = ADRIA.DataCube(s_rac[:, NORTH_MASK, :]; timesteps=ref_years, sites=1:count(NORTH_MASK), scenarios=1:1)
central_res = ADRIA.DataCube(s_rac[:, CENTRAL_MASK, :]; timesteps=ref_years, sites=1:count(CENTRAL_MASK), scenarios=1:1)
south_res = ADRIA.DataCube(s_rac[:, SOUTH_MASK, :]; timesteps=ref_years, sites=1:count(SOUTH_MASK), scenarios=1:1)


f = Figure(; Dict{Symbol,Any}(:size => (1600, 1600))...)
ax1 = plot_region(
    f,
    1,
    1,
    "North GBR",
    ltmp_north,
    north_res
)
ax2 = plot_region(
    f,
    1,
    2,
    "Central GBR",
    ltmp_central,
    central_res
)
ax3 = plot_region(
    f,
    2,
    1,
    "South GBR",
    ltmp_south,
    south_res;
    showlegend=true,
    legend_row=2,
    legend_col=2
)
linkyaxes!(ax1, ax2, ax3)

resize_to_layout!(f)
