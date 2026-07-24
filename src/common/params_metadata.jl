"""
Sidecar `.meta.toml` metadata for serialised calibration-result parameter vectors.

Records enough to detect and transparently upgrade parameter vectors serialised under an
older `CalibConfig` schema (see `PARAM_SCHEMA_VERSION`), so old `results.dat`/
`intermediate_results.dat` files keep working after the parameter space grows.
"""

function _meta_path(dat_path::String)::String
    return replace(dat_path, r"\.dat$" => ".meta.toml")
end

"""
    write_params_metadata(dat_path::String, n_params::Int)::Nothing

Write a `.meta.toml` sidecar next to `dat_path` recording the parameter count, current
`PARAM_SCHEMA_VERSION`, and the `CoralBloxCalib` package version that produced it, so
`load_calibrated_params` can later tell whether a serialised result predates a
parameter-space change.
"""
function write_params_metadata(dat_path::String, n_params::Int)::Nothing
    meta = Dict(
        "n_params" => n_params,
        "param_schema_version" => PARAM_SCHEMA_VERSION,
        "package_version" => string(pkgversion(@__MODULE__)),
        "created_at" => string(now(UTC)),
    )
    open(_meta_path(dat_path), "w") do io
        TOML.print(io, meta)
    end
    return nothing
end

"""
    load_calibrated_params(dat_path::String, expected_n::Int)::Vector{Float64}

Load a serialised calibration-result parameter vector, transparently upgrading vectors
shorter than `expected_n` by padding missing trailing values with `1.0`.

`1.0` is only a safe pad value because the current (and so far only) schema growth is the
trailing `dist_std` group-scale block, for which `1.0` is a no-op (see `DIST_STD_SCALE_LB`).
A future parameter-space change whose no-op value isn't `1.0` will need this reconsidered.

Reads the `.meta.toml` sidecar (if present, written by `write_params_metadata`) purely to
report which schema version produced a too-short vector; falls back to a raw length
comparison (with a less specific warning) for results saved before this sidecar existed.
"""
function load_calibrated_params(dat_path::String, expected_n::Int)::Vector{Float64}
    params = deserialize(dat_path)
    n = length(params)

    n == expected_n && return params

    if n > expected_n
        throw(ArgumentError(
            "$dat_path has $n params, more than the current schema's $expected_n - " *
            "can't safely truncate."
        ))
    end

    meta_path = _meta_path(dat_path)
    version_str = if isfile(meta_path)
        meta = TOML.parsefile(meta_path)
        "param_schema_version=$(get(meta, "param_schema_version", "unknown")), " *
        "package_version=$(get(meta, "package_version", "unknown")) (from $meta_path)"
    else
        "no .meta.toml sidecar found"
    end
    @warn "Loaded $n-param result ($version_str), padding $(expected_n - n) missing " *
          "trailing value(s) with 1.0 (no-op scale) to match the current $expected_n-param " *
          "schema. Path: $dat_path"

    return vcat(params, fill(1.0, expected_n - n))
end
