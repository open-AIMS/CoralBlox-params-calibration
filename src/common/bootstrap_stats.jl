using Bootstrap
using Random

"""
    BlockSampling(nrun)

Circular (moving) block resampling for [`block_bootstrap_ci`](@ref). `Bootstrap.jl`
itself only ships iid-style resamplers (`BasicSampling`, `BalancedSampling`, ...) plus a
maximum-entropy resampler for dependent series that regenerates synthetic values from a
single series' order statistics rather than resampling paired rows — neither respects the
year-to-year block structure we need here, so this subtypes `Bootstrap.jl`'s (unexported)
`BootstrapSampling` to plug a custom resampling scheme into its `NonParametricBootstrapSample`/
`confint`/`stderror` machinery.
"""
struct BlockSampling <: Bootstrap.BootstrapSampling
    nrun::Int
end

function _block_indices(n::Int, block_size_range::UnitRange{Int})::Vector{Int}
    idxs = Int[]
    while length(idxs) < n
        block_size = rand(block_size_range)
        start = rand(1:n)
        for k in 0:(block_size - 1)
            length(idxs) == n && break
            push!(idxs, mod(start + k - 1, n) + 1)
        end
    end
    return idxs
end

_iid_indices(n::Int)::Vector{Int} = rand(1:n, n)

"""
    _bootstrap_replicates(indices_fn, statistic, data, B)::Vector{Float64}

Draw `B` bootstrap replicates of `statistic(rows)`, rejecting and redrawing any resample
for which `statistic` isn't finite. Needed for ratio-type statistics (e.g. Nash-Sutcliffe
efficiency) on small samples: a resample can legitimately draw the same row `n` times,
making the resampled data's own variance exactly zero — a 0/0 or x/0 case the statistic
itself can't meaningfully resolve. Injecting a sentinel value for these would distort the
CI/SE; excluding them (a handful of retries, not a materially smaller effective `B`) is
the standard treatment for degenerate bootstrap resamples of this kind.
"""
function _bootstrap_replicates(
    indices_fn::Function, statistic::Function, data::AbstractMatrix, B::Int
)::Vector{Float64}
    n = size(data, 1)
    t1 = Vector{Float64}(undef, B)
    b = 1
    while b <= B
        t = statistic(@view data[indices_fn(n), :])
        isfinite(t) || continue
        t1[b] = t
        b += 1
    end
    return t1
end

"""
    block_bootstrap_ci(statistic, data::AbstractMatrix; B=2000, block_size_range=2:3,
                        ci_level=0.95)::NamedTuple

Circular block bootstrap CI for `statistic(rows)`, resampling contiguous (wrap-around)
row-blocks of `data` to respect row-order dependence (e.g. a time axis) — rows stay
intact, so paired columns (e.g. matched sim/obs values) resample together. Builds a
`Bootstrap.NonParametricBootstrapSample` by hand (neither the block resampling nor the
degenerate-resample rejection in [`_bootstrap_replicates`](@ref) is provided by
`Bootstrap.jl`), then defers to the package's own `confint`/`stderror` for the
interval/SE computation. Returns a `NamedTuple` with `estimate`, `lo`, `hi`, `se`.
"""
function block_bootstrap_ci(
    statistic::Function,
    data::AbstractMatrix;
    B::Int=2000,
    block_size_range::UnitRange{Int}=2:3,
    ci_level::Float64=0.95,
)::NamedTuple
    t0 = statistic(data)
    t1 = _bootstrap_replicates(n -> _block_indices(n, block_size_range), statistic, data, B)

    bs = NonParametricBootstrapSample((t0,), (t1,), statistic, data, BlockSampling(B))
    _, lo, hi = only(confint(bs, PercentileConfInt(ci_level)))
    return (estimate=t0, lo=lo, hi=hi, se=only(stderror(bs)))
end

"""
    iid_bootstrap_ci(statistic, data::AbstractMatrix; B=2000, ci_level=0.95)::NamedTuple

Iid bootstrap CI for `statistic(rows)`, resampling rows of `data` with replacement
(preserves row-pairing, e.g. matched sim/obs columns). Uses `Bootstrap.jl`'s own
`NonParametricBootstrapSample`/`confint`/`stderror` for the interval/SE computation, but
draws replicates via [`_bootstrap_replicates`](@ref) rather than `Bootstrap.bootstrap`
directly — its built-in driver has no hook to reject degenerate resamples (see
[`_bootstrap_replicates`](@ref)), which matters here since the reefs this is used for
have as few as 4 observed years, making an all-identical resample non-negligibly likely.
Returns a `NamedTuple` with `estimate`, `lo`, `hi`, `se`.
"""
function iid_bootstrap_ci(
    statistic::Function,
    data::AbstractMatrix;
    B::Int=2000,
    ci_level::Float64=0.95,
)::NamedTuple
    t0 = statistic(data)
    t1 = _bootstrap_replicates(_iid_indices, statistic, data, B)

    bs = NonParametricBootstrapSample((t0,), (t1,), statistic, data, BasicSampling(B))
    _, lo, hi = only(confint(bs, PercentileConfInt(ci_level)))
    return (estimate=t0, lo=lo, hi=hi, se=only(stderror(bs)))
end

"""
    bootstrap_median_ci(point_estimates; B=2000, ci_level=0.95)::NamedTuple

Second-level bootstrap over a group of point estimates: resample which estimates are
included (with replacement, same count as `point_estimates`), take the median of each
resample, repeat `B` times. Thin wrapper around `Bootstrap.bootstrap(median, ...)` — the
inputs here are always finite per-reef point estimates, so no rejection sampling is
needed. Returns a `NamedTuple` with `median`, `lo`, `hi`, `se`.
"""
function bootstrap_median_ci(
    point_estimates::AbstractVector{<:Real};
    B::Int=2000,
    ci_level::Float64=0.95,
)::NamedTuple
    bs = bootstrap(median, point_estimates, BasicSampling(B))
    _, lo, hi = only(confint(bs, PercentileConfInt(ci_level)))
    return (median=only(original(bs)), lo=lo, hi=hi, se=only(stderror(bs)))
end

"""
    bootstrap_mean_ci(point_estimates; B=2000, ci_level=0.95)::NamedTuple

Second-level bootstrap over a group of point estimates: resample which estimates are
included (with replacement, same count as `point_estimates`), take the mean of each
resample, repeat `B` times. Same resampling scheme as [`bootstrap_median_ci`](@ref), with
`mean` in place of `median` — used for aggregates where the magnitude (not the signed
central tendency) of the point estimates is of interest, e.g. mean absolute bias. Returns
a `NamedTuple` with `mean`, `lo`, `hi`, `se`.
"""
function bootstrap_mean_ci(
    point_estimates::AbstractVector{<:Real};
    B::Int=2000,
    ci_level::Float64=0.95,
)::NamedTuple
    bs = bootstrap(mean, point_estimates, BasicSampling(B))
    _, lo, hi = only(confint(bs, PercentileConfInt(ci_level)))
    return (mean=only(original(bs)), lo=lo, hi=hi, se=only(stderror(bs)))
end
