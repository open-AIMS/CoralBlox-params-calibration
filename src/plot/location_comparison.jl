"""
    plot_location_comparison(dom, raw_data::Array{Float64,3}, reef_id::String, cyc_scens, dhw_scens, disturbances; observations::LocationDataStore=COMBINED_STORE)::Figure
    plot_location_comparison(dom, raw_data::Array{Float64,3}, ltmp_loc_index::Int64, cyc_scens, dhw_scens, disturbances; observations::LocationDataStore=COMBINED_STORE)::Figure

Plot LTMP Manta Tow Coral Cover against the modelled cover output given the LTMP location index.
"""
function plot_location_comparison(
    raw_data::Array{Float64,3},
    reef_id::String,
    dhw_scens::AbstractMatrix{T},
    cyc_scens::AbstractArray{Float64,3},
    disturbances;
    observations::LocationDataStore=COMBINED_STORE,
    hide_disturbances=false
)::Figure where {T<:Real}
    ltmp_loc_index = findfirst(VALIDATION_STORE.ltmp_unique_ids .== reef_id)
    return plot_location_comparison(
        raw_data, ltmp_loc_index, cyc_scens, dhw_scens, disturbances;
        observations=observations
    )
end
function plot_location_comparison(
    raw_data::Array{Float64,3},
    ltmp_loc_idx::Int64,
    dhw_scens::AbstractMatrix{T},
    cyc_scens::AbstractArray{Float64,3},
    disturbances;
    observations::LocationDataStore=COMBINED_STORE,
    hide_model_disturbances=false,
    hide_ltmp_disturbances=false,
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Figure where {T<:Real}
    f = Figure(; size=(1000, 1000))

    title_text = _location_err_title(
        raw_data, ltmp_loc_idx;
        fig_opts,
        observations=observations
    )

    # Plot LTMP and CoralBlox covers
    # cover_limits = (xlimits, (0, maximum(raw_data) * 1.1))
    current_row = 1
    modelled_ltmp_axis_opts = Dict(
        # :title => title_text,
        :titlesize => FONT_SIZES.title,
        :ylabel => "Relative coral cover",
        BASE_AXIS_OPTS...
    )

    plot_modelled_v_ltmp!(
        f, current_row, raw_data, ltmp_loc_idx;
        observations=observations, axis_opts=modelled_ltmp_axis_opts
    )
    legend_modelled_v_ltmp!(f, current_row)
    current_row += 1

    # These are RME_UNIQUE_ID. Need to select locations from LocationDataStore
    reef_id = observations.ltmp_unique_ids[ltmp_loc_idx]

    # Need to use RME_GBRMPA_ID to select locations from dhw_scens and cyc_scens
    spatial_data = observations.domain_gpkg[:, [:RME_UNIQUE_ID, :RME_GBRMPA_ID]]
    reef_gbrmpa_id = spatial_data[spatial_data.RME_UNIQUE_ID.==reef_id, :RME_GBRMPA_ID][1]

    if !hide_model_disturbances
        ax_opts_coralblox_disturbances = Dict(
            BASE_AXIS_OPTS...,
            :ylabel => "DHW",
            :height => 100
        )
        loc_dhw_scens = dhw_scens[locs=At(reef_gbrmpa_id)]
        loc_cyc_scens = cyc_scens[locs=At(reef_gbrmpa_id)]

        plot_coralblox_disturbances!(
            f, current_row, loc_dhw_scens, loc_cyc_scens;
            axis_opts=ax_opts_coralblox_disturbances,
        )
        legend_coralblox_disturbances!(f, current_row)
        current_row += 1
    end

    # Plot CoralBlox Observed DHW and Cyclone Categories
    if !hide_ltmp_disturbances
        if reef_id ∈ disturbances.locations
            loc_ltmp_disturbances = disturbances[locations=At(reef_id)]
            plot_ltmp_disturbances!(
                f, current_row, loc_ltmp_disturbances;
                axis_opts=ax_opts_coralblox_disturbances
            )
            ax_row_taxa = 4

            legend_ltmp_disturbances!(f, current_row, loc_ltmp_disturbances)
            current_row += 1
        else
            @warn("Reef $reef_id / $ltmp_loc_idx not found in LTMP disturbances DataFrame")
            ax_row_taxa = 3
        end
    end

    cb_loc_id = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    cover = rs_raw.raw[:, :, cb_loc_id]
    plot_taxa_props!(f, current_row, cover; ax_opts=BASE_AXIS_OPTS)

    legend_taxa_props!(f, current_row)

    titlesize = get(fig_opts, :titlesize, 20)
    Label(f[0, :], title_text, fontsize=titlesize)

    return f
end

# TODO Rename this. Maybe plot_location_comparison_highlights or something
function plot_location_comparison_highlights(
    raw_data::Array{Float64,3},
    reef_ids::Vector{String},
    cyc_scens,
    dhw_scens,
    disturbances;
    observations::LocationDataStore=COMBINED_STORE
)::Figure
    # Extra information specific to the given location
    # reef_id = observations.ltmp_unique_ids[ltmp_loc_idx]

    f = Figure(; size=(1500, 1500), title="Model prediction versus historic data")
    ax_11 = f[1, 1] = GridLayout()
    ax_12 = f[1, 2] = GridLayout()
    ax_21 = f[2, 1] = GridLayout()
    ax_22 = f[2, 2] = GridLayout()

    # TODO try legend on the bottom/horizontal
    ax_legend = f[1:end, 3] = GridLayout()

    ltmp_dist_colors = getproperty(Makie.ColorSchemes, COLORMAPS.ltmp_disturbances)
    ltmp_dist_labels = collect(disturbances.disturbances)
    taxa_colors = getproperty(Makie.ColorSchemes, COLORMAPS.taxa)
    taxa = String.(ADRIA.functional_group_names())
    taxa_labels = titlecase.(join.(split.(taxa, "_"), " "))

    # Legend: LTMP and CoralBlox
    ltmp_el = MarkerElement(color=:red, marker=:circle, markersize=10)
    model_el = LineElement(color=:red, linewidth=2)

    # Legend: Cyc and DHW
    dhw_el = LineElement(color=:orange, linewidth=4, linestyle=:dash)
    cyc_ele = LineElement(color=:black, linewidth=4, linestyle=:dash)
    dist_labels = ["DHW", "Cyclone / COTS"]

    # Legend: LTMP Disturbances
    ltmp_dist_ele = [
        LineElement(color=ltmp_dist_colors[i], linewidth=2)
        for i in 1:length(ltmp_dist_colors)
    ]

    # Legend: Functional groups
    taxa_ele = [LineElement(color=taxa_colors[i], linewidth=2) for i in 1:length(taxa_colors)]

    Legend(
        ax_legend[1:end, 1],
        [[ltmp_el, model_el], vcat(dhw_el, cyc_ele), ltmp_dist_ele, taxa_ele],
        [["LTMP", "CoralBlox"], dist_labels, ltmp_dist_labels, taxa_labels],
        ["Coral Cover", "CoralBlox Disturbances", "LTMP Disturbances", "Functional groups"],
        framevisible=false,
        valign=:top,
    )

    axis = [ax_11, ax_12, ax_21, ax_22]

    axis_title = ["(a)", "(b)", "(c)", "(d)"]

    for (ax_idx, ax) in enumerate(axis)
        reef_id = reef_ids[ax_idx]

        loc_idx = findfirst(observations.ltmp_unique_ids .== reef_id)
        title_text = _location_err_title(
            raw_data,
            loc_idx;
            observations=observations
        )

        # Plot LTMP and CoralBlox covers
        # cover_limits = (xlimits, (0, maximum(raw_data) * 1.1))
        modelled_v_ltmp_row = 1
        modelled_ltmp_axis_opts = Dict(
            BASE_AXIS_OPTS...,
            :title => axis_title[ax_idx],
            :ylabel => "Relative total cover",
            :titlesize => FONT_SIZES.title,
        )
        obs_modelled_v_ltmp, sim_modelled_v_ltmp = plot_modelled_v_ltmp!(
            ax, modelled_v_ltmp_row, raw_data, loc_idx;
            observations=observations, axis_opts=modelled_ltmp_axis_opts
        )

        # Plot CoralBlox Observed DHW and Cyclone Categories
        ax_opts_coralblox_disturbances = Dict(
            BASE_AXIS_OPTS...,
            :ylabel => "DHW",
            :height => 100
        )
        ax_row_coralblox_disturbances::Int64 = 2
        loc_dhw_scens::YAXArray{Float32,1} = dhw_scens[locs=At(reef_id)]
        loc_cyc_scens::YAXArray{Float64,2} = cyc_scens[locations=At(reef_id)]
        plot_coralblox_disturbances!(
            ax, ax_row_coralblox_disturbances, loc_dhw_scens, loc_cyc_scens;
            axis_opts=ax_opts_coralblox_disturbances,
        )

        if reef_id ∈ disturbances.locations
            loc_ltmp_disturbances = disturbances[locations=At(reef_id)]
            ax_row_ltmp_disturbances = 3
            ax_opts_ltmp_disturbances = Dict(BASE_AXIS_OPTS..., :height => 100)
            plot_ltmp_disturbances!(
                ax, ax_row_ltmp_disturbances, loc_ltmp_disturbances;
                axis_opts=ax_opts_ltmp_disturbances
            )
            ax_row_taxa = 4
        else
            @warn("Reef $reef_id / $ltmp_loc_idx not found in LTMP disturbances DataFrame")
            ax_row_taxa = 3
        end

        ltmp_loc_idx = findfirst(VALIDATION_STORE.ltmp_unique_ids .== reef_id)
        cb_loc_id = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
        cover = rs_raw.raw[:, :, cb_loc_id]
        plot_taxa_props!(ax, ax_row_taxa, cover; ax_opts=BASE_AXIS_OPTS)
    end

    return f
end

function plot_modelled_v_ltmp!(
    fig::Union{Figure,GridLayout},
    ax_row::Int64,
    raw_data,
    ltmp_loc_idx;
    observations::LocationDataStore=COMBINED_STORE,
    dom=dom,
    axis_opts::Dict=Dict()
)
    ax = Axis(fig[ax_row, 1]; axis_opts...)
    return plot_modelled_v_ltmp!(ax, raw_data, ltmp_loc_idx; observations=observations)
end
function plot_modelled_v_ltmp!(
    ax::Axis,
    raw_data,
    ltmp_loc_idx;
    observations::LocationDataStore=COMBINED_STORE,
    dom=dom
)
    domain_loc_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    loc_k_area = ADRIA.site_k_area(dom)[domain_loc_idx]
    loc_area = ADRIA.loc_area(dom)[domain_loc_idx]

    loc_cover = dropdims(
        sum(raw_data[:, :, domain_loc_idx], dims=2),
        dims=2) .* loc_k_area ./ loc_area

    obs_loc_data = observations.ltmp_coral_cover[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))
    obs_tf = (START_YEAR:END_YEAR)[not_missing_obs]

    obs = scatter!(ax, obs_tf, obs_loc_data[not_missing_obs], color=(:red, 0.5), markersize=20)
    sim = lines!(ax, 2008:2022, loc_cover, color=:red)
    return obs, sim
end

function legend_modelled_v_ltmp!(f::Figure, row::Int64; col::Int64=2)
    ltmp_el = MarkerElement(color=:red, marker=:circle, markersize=10)
    model_el = LineElement(color=:red, linewidth=2)
    Legend(
        f[row, col],
        [ltmp_el, model_el],
        ["LTMP", "CoralBlox"],
        "Coral Cover",
        valign=:top,
    )
end

function plot_coralblox_disturbances!(
    fig::Union{Figure,GridLayout},
    ax_row::Int64,
    loc_dhw_scens::AbstractVector{R1},
    loc_cyclone_scens::AbstractMatrix{R2};
    axis_opts::Dict=Dict()
)::Nothing where {R1,R2<:Real}
    ax_dhw = Axis(fig[ax_row, 1]; axis_opts...)
    plot_dhw_scens!(ax_dhw, loc_dhw_scens)

    ax_cyclone = Axis(fig[ax_row, 1]; axis_opts...)
    hidespines!(ax_cyclone)
    hidedecorations!(ax_cyclone)
    plot_cyclone_scens!(ax_cyclone, loc_cyclone_scens)

    return nothing
end

function legend_coralblox_disturbances!(f::Figure, row::Int64; col::Int64=2)
    dhw_el = LineElement(linewidth=4, color=:orange, linestyle=:dot)
    dist_el = LineElement(linewidth=4, color=:black, linestyle=:dash)
    return Legend(
        f[row, col], [dhw_el, dist_el], ["DHW", "Cyclone/COTS"], "Disturbances", valign=:top
    )
end

function plot_dhw_scens!(ax::Axis, dhw_scens::AbstractVector{R}) where {R<:Real}
    dhw_data::Vector{Float64} = Float64.(collect(dhw_scens))
    timesteps::Vector{Int64} = collect(dhw_scens.timesteps.val.data)
    return lines!(ax, timesteps, dhw_data, color=:orange, linestyle=:dot, linewidth=4)
end

function plot_cyclone_scens!(ax::Axis, cyc_scens::AbstractMatrix{Float64})::Nothing
    sum_cyc_scens = dropdims(sum(cyc_scens; dims=2), dims=2)
    target_years = sum_cyc_scens[sum_cyc_scens.>0].timesteps.val.data
    vlines!.(ax, target_years; ymin=0, color=:black, linestyle=:dash, linewidth=3)
    return nothing
end

function legend_ltmp_disturbances!(
    f::Figure, row::Int64, loc_disturbances::YAXArray{Int16,2}; col::Int64=2
)
    disturbance_mask = read(dropdims(sum(loc_disturbances, dims=1), dims=1) .> 0)
    colors = Makie.ColorSchemes.seaborn_bright6[disturbance_mask]
    labels = collect(loc_disturbances.disturbances.val)[disturbance_mask]
    return legend_ltmp_disturbances!(f, row, colors, labels; col=col)
end
function legend_ltmp_disturbances!(
    f::Figure, row::Int64, disturbances::YAXArray{Int16,3}; col::Int64=2
)
    colors = Makie.ColorSchemes.seaborn_bright6
    labels = collect(disturbances.disturbance)
    return legend_ltmp_disturbances!(f, row, colors, labels; col=col)
end
function legend_ltmp_disturbances!(
    f::Figure, row::Int64, colors, labels::Vector{String};
    col::Int64=2
)
    n_colors = length(colors)
    ltmp_dist_ele = [LineElement(color=colors[i], linewidth=2) for i in 1:n_colors]
    return Legend(f[row, col], ltmp_dist_ele, labels, "LTMP Disturbances", valign=:top)
end

function legend_taxa_props!(f::Figure, row::Int64; col=2)
    taxa_colors = Makie.ColorSchemes.seaborn_colorblind6[1:5]
    n_taxa = length(taxa_colors)
    taxa = String.(ADRIA.functional_group_names())
    taxa_labels = titlecase.(join.(split.(taxa, "_"), " "))
    taxa_ele = [LineElement(color=taxa_colors[i], linewidth=2) for i in 1:n_taxa]

    return Legend(f[row, col], taxa_ele, taxa_labels, "Functional groups", valign=:top,)
end

"""
    plot_observation_locs(
        calibration_store::LocationDataStore, validation_store::LocationDataStore
    )::Figure

Map showing calibration and valiration locations with different colors.
"""
function plot_observation_locs(
    calibration_store::LocationDataStore, validation_store::LocationDataStore
)::Figure
    domain_gpkg = calibration_store.domain_gpkg
    calib_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈Ref(calibration_store.ltmp_unique_ids), :]
    valid_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈Ref(validation_store.ltmp_unique_ids), :]

    f = Figure(size=FIG_SIZE[:map])
    ax = GeoAxis(
        f[1, 1],
        dest="+proj=latlong +datum=WGS84",
        title="Observation Locations",
        titlesize=FONT_SIZES[:title]
    )
    try
        poly!(ax, domain_gpkg.geom, color=:black)
    catch
        poly!(ax, domain_gpkg.geometry, color=:black)
    end

    obs = (Calibration=:red, Validation=:blue)

    try
        scatter!(ax, calib_gpkg.X_COORD, calib_gpkg.Y_COORD; markersize=10, color=obs.Calibration, alpha=0.5)
        scatter!(ax, valid_gpkg.X_COORD, valid_gpkg.Y_COORD; markersize=10, color=obs.Validation, alpha=0.5)
    catch
        scatter!(ax, calib_gpkg.LON, calib_gpkg.LAT; markersize=10, color=obs.Calibration, alpha=0.5)
        scatter!(ax, valid_gpkg.LON, valid_gpkg.LAT; markersize=10, color=obs.Validation, alpha=0.5)
    end


    els = [MarkerElement(color=c, marker=:circle, markersize=10) for c in values(obs)]
    labels = collect(string.(keys(obs)))
    @info(els)
    @info(labels)
    Legend(
        f[2, 1],
        els,
        labels,
        valign=:top,
        orientation=:horizontal
    )

    return f
end
