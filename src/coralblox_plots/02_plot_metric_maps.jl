metric_map_fig = begin
    fig = Figure(; size=(550, 650), fontsize=9pt)

    markersize = 10
    strokecolor = (:gray, 0.5)
    grid_color = RGBAf(0, 0, 0, 0.10)
    rmse_colormap = :roma
    srcc_colormap = :cool

    srcc_colorrange = (-1, 1)
    srcc_ticks_range = srcc_colorrange[1]:0.5:srcc_colorrange[2]
    srcc_colorbar_ticks = (srcc_ticks_range, string.(srcc_ticks_range))

    rmse_colorrange = (-0.2, 0.2)
    rmse_ticks_range = rmse_colorrange[1]:0.1:rmse_colorrange[2]
    rmse_colorbar_ticks = (rmse_ticks_range, string.(rmse_ticks_range))

    g0 = fig[1, 1:2] = GridLayout()
    g1 = fig[2, 1:2] = GridLayout()
    g11 = g1[1, 1] = GridLayout()
    g12 = g1[1, 2] = GridLayout()
    gr2 = fig[3, 1:2] = GridLayout()
    g21 = gr2[1, 1] = GridLayout()
    g22 = gr2[1, 2] = GridLayout()

    plot_metric_map!(
        g11,
        c_rmse_diffs,
        geometries,
        c_lon,
        c_lat;
        axis_opts=Dict(
            :title => "",
            :titlesize => 9pt,
            :xlabelsize => 9pt,
            :ylabelsize => 9pt,
            :yticklabelsize => 9pt,
            :xticklabelsize => 9pt,
            :xgridcolor => grid_color,
            :ygridcolor => grid_color,
            :xticklabelsvisible => false),
        opts=Dict(
            :colorrange => (-0.2, 0.2),
            :colorbar_label => "",
            :colorbar_ticklabelsize => 9pt,
            :colorbar_ticks => rmse_colorbar_ticks,
            :markersize => markersize,
            :alpha => 0.8,
            :strokewidth => 0.1,
            :strokecolor => strokecolor,
            :colormap => rmse_colormap,
            :colorbar_visible => false,
        )
    )

    plot_metric_map!(
        g12,
        collect(c_srcc_),
        geometries,
        c_lon,
        c_lat;
        axis_opts=Dict(
            :title => "",
            :titlesize => 9pt,
            :xlabelsize => 9pt,
            :ylabelsize => 9pt,
            :yticklabelsize => 9pt,
            :xticklabelsize => 9pt,
            :xgridcolor => grid_color,
            :ygridcolor => grid_color,
            :xticklabelsvisible => false,
            :yticklabelsvisible => false
        ),
        opts=Dict(
            :colorbar_label => "",
            :colorbar_ticklabelsize => 9pt,
            :colorbar_ticks => srcc_colorbar_ticks,
            :markersize => markersize,
            :alpha => 0.6,
            :strokewidth => 0.1,
            :strokecolor => strokecolor,
            :colormap => srcc_colormap,
            :colorbar_visible => false,
        )
    )

    plot_metric_map!(
        g21,
        v_rmse_diffs,
        geometries,
        v_lon,
        v_lat;
        axis_opts=Dict(
            :title => "",
            :titlesize => 9pt,
            :xlabelsize => 9pt,
            :ylabelsize => 9pt,
            :yticklabelsize => 9pt,
            :xticklabelsize => 9pt,
            :xgridcolor => grid_color,
            :ygridcolor => grid_color,
            :valign => :top,
        ),
        opts=Dict(
            :colorrange => rmse_colorrange,
            :colorbar_label => "",
            :colorbar_ticklabelsize => 9pt,
            :colorbar_ticks => rmse_colorbar_ticks,
            :markersize => markersize,
            :alpha => 0.8,
            :strokewidth => 0.1,
            :strokecolor => strokecolor,
            :colormap => rmse_colormap,
            :colorbar_visible => false,
        )
    )

    plot_metric_map!(
        g22,
        collect(v_srcc_),
        geometries,
        v_lon,
        v_lat;
        axis_opts=Dict(
            :title => "",
            :titlesize => 9pt,
            :xlabelsize => 9pt,
            :ylabelsize => 9pt,
            :yticklabelsize => 9pt,
            :xticklabelsize => 9pt,
            :xgridcolor => grid_color,
            :ygridcolor => grid_color,
            :valign => :top,
            :yticklabelsvisible => false
        ),
        opts=Dict(
            :colormap => srcc_colormap,
            :colorbar_label => "",
            :colorbar_ticklabelsize => 9pt,
            :colorbar_ticks => srcc_colorbar_ticks,
            :markersize => markersize,
            :alpha => 0.6,
            :strokewidth => 0.1,
            :strokecolor => strokecolor,
            :colorbar_visible => false,
        )
    )

    lon_padding = (0, 0, 0, 25)
    lat_padding = (0, -20, 0, 0)
    left_label_padding = (0, 20, 0, 0)
    top_label_padding = (0, 0, 60, -30)

    gcb1 = g0[1, 1] = GridLayout()
    gcb2 = g0[1, 2] = GridLayout()
    cb_rmse = Colorbar(
        gcb1[1, 1];
        colorrange=rmse_colorrange,
        colormap=rmse_colormap,
        ticks=rmse_colorbar_ticks,
        ticklabelsize=9pt,
    )
    Label(gcb1[1, 1], "ΔRMSE (Benchmark - Model) ", padding=(0, 0, 60, -30), fontsize=11pt, justification=:center, halign=:center)

    cb_srcc = Colorbar(
        gcb2[1, 1];
        colorrange=srcc_colorrange,
        colormap=srcc_colormap,
        ticks=srcc_colorbar_ticks,
        ticklabelsize=9pt,)
    Label(gcb2[1, 1], "SRCC", padding=(0, 0, 60, -30), fontsize=11pt, justification=:center, halign=:center)

    Label(g21[1, 1, Bottom()], "Longitude", padding=lon_padding, fontsize=9pt)
    Label(g22[1, 1, Bottom()], "Longitude", padding=lon_padding, fontsize=9pt)

    Label(g11[1, 1, Left()], "Latitude", rotation=π / 2, padding=lat_padding, fontsize=9pt, tellwidth=false)
    Label(g21[1, 1, Left()], "Latitude", rotation=π / 2, padding=lat_padding, fontsize=9pt, tellwidth=false)

    Label(g11[1, 1, Left()], "Calibration reefs", rotation=π / 2, padding=left_label_padding, fontsize=11pt,)
    Label(g21[1, 1, Left()], "Validation reefs", rotation=π / 2, padding=left_label_padding, fontsize=11pt,)

    colgap!(g0, 32)
    colgap!(g1, -60)
    colgap!(gr2, -60)
    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, -10)

    cb_srcc.vertical, cb_rmse.vertical = false, false
    cb_srcc.labelvisible, cb_rmse.labelvisible = false, false
    cb_srcc.width, cb_rmse.width = 180, 180

    resize_to_layout!(fig)
    fig
end
save(fig_path * "/metric_maps.png", metric_map_fig; px_per_unit=(300 / inch))