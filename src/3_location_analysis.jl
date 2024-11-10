if !isdefined(Main, :reload_calibration_results)
    reload_calibration_results = true
end

if !isdefined(Main, :results_setup) || reload_calibration_results
    include("results_analysis.jl")
    results_setup = true
    reload_calibration_results = false
end

function plot_coral_param(
    loc::String,
    param_name::String,
    category,
    group,
    data
)::Figure

    fig = Figure(; size=(1300, 900))
    ax = Axis(
        fig[1,1],
        xticks = 1:length(category) / 5,
        title = "$(location_unique_id): $(param_name)"
    )

    barplot!(
        ax,
        category,
        data,
        dodge = group,
        color = Makie.wong_colors()[group]
    )

    labels = String.(ADRIA.functional_group_names())
    elements = [PolyElement(polycolor = Makie.wong_colors()[i]) for i in 1:length(labels)]
    title = "Functional Groups"

    Legend(fig[1,2], elements, labels, title)

    return fig
end

function plot_taxa_props(
    loc::String,
    cover;
)::Figure
    cover = reshape(cover, (15, 7, 5))
    cover = dropdims(sum(cover, dims=2), dims=2) ./ dropdims(sum(cover, dims=(2, 3)), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022
    f = Figure(; size=(1300, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Proportion",
        title="$(loc): Functional Group Cover Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, cover, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)

    return f
end

function plot_bleaching_mortality(
    loc::String,
    loc_cover,
    loc_bleaching
)::Figure
    prop_sc_fg_cover = permutedims(reshape(loc_cover, (15, 7, 5)), (1, 3, 2))
    prop_sc_fg_cover = prop_sc_fg_cover ./ sum(prop_sc_fg_cover, dims=(2, 3))
    perc_loss = dropdims(sum(prop_sc_fg_cover .* loc_bleaching, dims=(2, 3)), dims=(2, 3))

    xs = 2008:2022
    f = Figure(; size=(1300, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Loss",
        title="$(loc): Proportion Cover Loss from Bleaching"
    )
    sr = lines!(xs, perc_loss)

    return f
end

function plot_cyclone_mortality(
    loc::String,
    loc_cover,
    loc_cyclone
)::Figure
    prop_sc_fg_cover = permutedims(reshape(loc_cover, (15, 7, 5)), (1, 3, 2))
    prop_sc_fg_cover = dropdims(sum(prop_sc_fg_cover ./ sum(prop_sc_fg_cover, dims=(2, 3)), dims=3), dims=3)
    perc_loss = dropdims(sum(prop_sc_fg_cover .* loc_cyclone, dims=2), dims=2)

    xs = 2008:2022
    f = Figure(; size=(1300, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Loss",
        title="$(loc): Proportion Cover Loss from Cyclones"
    )
    sr = lines!(xs, perc_loss)
    return f
end

function plot_target_proportions(loc::String, reefmod_taxa)::Figure
    xs = 2008:2022
    normalised_comp = reefmod_taxa ./ sum(reefmod_taxa, dims = 2)
    non_missing_mask = (!).(ismissing.(reefmod_taxa[:, 1]))

    f = Figure(;size=(1300, 900))

    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Loss",
        xticks = xs,
        title="$(loc): Target Coral Composition",
        limits=(nothing, nothing, 0, 1)
    )

    xs_m = xs[non_missing_mask]
    dt = permutedims(Float64.(normalised_comp[non_missing_mask, :]), (2, 1))

    sr = series!(xs_m, dt, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)
    return f
end

function error_functions(
    raw_data,
    ltmp_loc_idx;
    obs_data=all_ltmp_reef,
    obs_idxs=all_ltmp_idxs,
    obs_loc_labels=ltmp_reef_data.RME_UNIQUE_ID,
    loc_k_areas=ADRIA.site_k_area(dom),
    loc_areas=ADRIA.loc_area(dom)
)
    loc_cover = dropdims(sum(raw_data, dims=2), dims=2) .* loc_k_areas' ./ loc_areas'
    if obs_idxs[ltmp_loc_idx] == -1
        return Figure()
    end

    obs_loc_data = obs_data[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))
    obs_tf = (start_year:end_year)[not_missing_obs]

    sim_data = loc_cover[:, obs_idxs[ltmp_loc_idx]]
    reef_id = obs_loc_labels[ltmp_loc_idx]


    rmse_::Float64 = rmse(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    cc_::Float64 = cor(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    maee_::Float64 = MAEE(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    bias_::Float64 = bias(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    @info "Location $(reef_id)"
    @info "RMSE:       $(trunc(rmse_, digits=4))"
    @info "Pearsons R: $(trunc(cc_, digits=4))"
    @info "MAEE:       $(trunc(maee_, digits=4))"
    @info "Bias:       $(trunc(bias_, digits=4))"
    return rmse_, cc_, maee_, bias_
end

corals = ADRIA.to_coral_spec(scens[1, :])

coral_cover = rs_raw.raw .* reshape(site_k_area(dom), (1, 1, 3806))

# Identify temporal regin of interest
temporal_range = 2008:2014

# Location coral params
limited_loc_pos = 4 # index of target location in
location_unique_id = limited_locations[limited_loc_pos]
domain_loc_pos = findfirst(x->x==location_unique_id, dom.loc_data.UNIQUE_ID)
ltmp_loc_pos = findfirst(x->!ismissing(x) && x==location_unique_id, ltmp_reef_data.RME_UNIQUE_ID)

println("---- Plotting Location $(location_unique_id) ----")
println("Reef Name: $(dom.loc_data.LOC_NAME_L[domain_loc_pos])")
println("Calibration Location Index: $(limited_loc_pos)")
println("Ltmp Location Index: $(ltmp_loc_pos)")

linear_ext = permutedims(reshape(corals.linear_extension, (7, 5)), (2, 1))
linear_ext[:, 7] .= 0.0
linear_ext .*= scale_factors[:, limited_loc_pos, 1]

mb_rate = permutedims(reshape(corals.mb_rate, (7, 5)), (2, 1))
mb_rate .*= scale_factors[:, limited_loc_pos, 2]

# Bar plot groupings
cat = [floor((i - 1) / 5) + 1 for i in 1:35]
grp = [(i - 1) % 5 + 1 for i in 1:35]

flt_linear_ext = reshape(linear_ext, (35,))

save_dir = "Outputs/$(prefix)/Location_$(location_unique_id)"

mkpath(save_dir)

f = plot_coral_param(location_unique_id, "Linear Extension", cat[1:30], grp[1:30], flt_linear_ext[1:30])
save("$(save_dir)/linexts.png", f)

flt_mb_rate = reshape(mb_rate, (35,))

f = plot_coral_param(location_unique_id, "Background Mortality", cat, grp, flt_mb_rate)
save("$(save_dir)/mbrates.png", f)

f = plot_taxa_props(location_unique_id, rs_raw.raw[:, :, domain_loc_pos])
save("$(save_dir)/taxa_props.png", f)

f = plot_bleaching_mortality(location_unique_id, rs_raw.raw[:, :, domain_loc_pos], rs_raw.bleaching_mortality[:, :, :, domain_loc_pos])
save("$(save_dir)/bleaching_mortality.png", f)

f = plot_cyclone_mortality(location_unique_id, rs_raw.raw[:, :, 1293], dom.cyclone_mortality_scens[:, 1293, :, scens[1, :cyclone_mortality_scenario]].data)
save("$(save_dir)/cyclone_mortality.png", f)

f = plot_target_proportions(location_unique_id, rm_ltmp_taxa[:, :, limited_loc_pos])
save("$(save_dir)/target_composition.png", f)

xs = dom.loc_data.X_COORD[domain_loc_pos]
ys = dom.loc_data.Y_COORD[domain_loc_pos]

f = ADRIA.viz.map(dom)
scatter!([xs], [ys], color=(:black, 0.0), markersize = 20, strokecolor=:black, strokewidth=2, overdraw=true)
save("$(save_dir)/loc_map.png", f)

location_comparison(rs_raw.raw, ltmp_loc_pos, save_dir)
error_functions(rs_raw.raw, ltmp_loc_pos)

loc_cov = rs_raw.raw[:, :, domain_loc_pos:domain_loc_pos]

f = temporal_size_class_proportions(loc_cov)
save("$(save_dir)/size_class_proportions.png", f)
