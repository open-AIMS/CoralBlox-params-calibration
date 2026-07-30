"""
Sidecar `.meta.toml` metadata for serialised calibration-result parameter vectors.

Records enough to detect parameter vectors serialised under an older `CalibConfig` schema
(see `PARAM_SCHEMA_VERSION`), so stale `results.dat`/`intermediate_results.dat` files are
rejected rather than silently reinterpreted when the parameter space changes.
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

Load a serialised calibration-result parameter vector, erroring if it was not written under
the current `PARAM_SCHEMA_VERSION`.

Vectors from older schemas cannot be upgraded in place. Schema v3 moved `dist_std` from a
trailing block of 5 group-level scale factors into the coral parameter block as 35 absolute
per-group/size-class values, which shifts the meaning of every index after the coral block.
Padding or truncating to match lengths would silently reinterpret unrelated parameters.

A vector of the expected length whose `.meta.toml` sidecar is missing is accepted with a
warning: every schema so far has had a distinct parameter count, so length alone identifies
it. A sidecar recording a *different* version is always fatal.
"""
function load_calibrated_params(dat_path::String, expected_n::Int)::Vector{Float64}
    params = deserialize(dat_path)
    n = length(params)

    meta_path = _meta_path(dat_path)
    version = isfile(meta_path) ?
        get(TOML.parsefile(meta_path), "param_schema_version", nothing) : nothing

    if n == expected_n && (version == PARAM_SCHEMA_VERSION || isnothing(version))
        isnothing(version) && @warn "No .meta.toml sidecar for $dat_path; accepting it as " *
            "schema v$PARAM_SCHEMA_VERSION on its $n-param length alone."
        return params
    end

    version_str = isnothing(version) ?
        "no .meta.toml sidecar found (predates schema versioning)" :
        "param_schema_version=$version"
    throw(ArgumentError(
        "$dat_path holds $n params ($version_str) but the current schema " *
        "v$PARAM_SCHEMA_VERSION expects $expected_n. Parameter vectors cannot be migrated " *
        "across schema versions - re-run the calibration from scratch."
    ))
end
