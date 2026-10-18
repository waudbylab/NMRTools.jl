module MakieExt

using NMRTools
using NMRTools: NMRData, AbstractNMRData, NMRDimension, FrequencyDimension,
                NonFrequencyDimension
using NMRTools: data, dims, refdims, refdims_title, label, scale, units,
                reorder, ForwardOrdered
using MulticomplexNumbers: Multicomplex, realest
using Statistics: quantile
using Makie
using Colors

include("MakieExt/plotutils.jl")
include("MakieExt/plot_1d.jl")
include("MakieExt/plot_2d.jl")
include("MakieExt/plot_pseudo2d.jl")
include("MakieExt/plot_3d.jl")

# Convenience: forward `nmrplot!(fig_like, spec)` to the underlying axis so
# the typical workflow `fig = nmrplot(spec); nmrplot!(fig, spec/2)` just
# works without manual destructuring.
NMRTools.nmrplot!(fap::Makie.FigureAxisPlot, args...; kwargs...) =
    NMRTools.nmrplot!(fap.axis, args...; kwargs...)

NMRTools.nmrplot!(ap::Makie.AxisPlot, args...; kwargs...) =
    NMRTools.nmrplot!(ap.axis, args...; kwargs...)

function NMRTools.nmrplot!(fig::Makie.Figure, args...; kwargs...)
    ax = Makie.current_axis(fig)
    isnothing(ax) && throw(ArgumentError("nmrplot!: Figure has no current axis. \
                                          Create one with `nmrplot(fig[1,1], spec)` \
                                          or pass an `Axis` directly."))
    return NMRTools.nmrplot!(ax, args...; kwargs...)
end

# `nmrplot!(spec; kwargs...)` with no axis argument falls back to Makie's
# `current_axis()` so the common workflow `nmrplot(spec); nmrplot!(spec/2)`
# works without an explicit handle.
function NMRTools.nmrplot!(spec::NMRData; kwargs...)
    ax = Makie.current_axis()
    isnothing(ax) && throw(ArgumentError("nmrplot!: no current axis. \
        Call nmrplot(spec) first, or pass an Axis / Axis3 explicitly."))
    return NMRTools.nmrplot!(ax, spec; kwargs...)
end

function NMRTools.nmrplot!(v::AbstractVector{<:NMRData}; kwargs...)
    ax = Makie.current_axis()
    isnothing(ax) && throw(ArgumentError("nmrplot!: no current axis. \
        Call nmrplot(specs) first, or pass an Axis / Axis3 explicitly."))
    return NMRTools.nmrplot!(ax, v; kwargs...)
end

# Fallback for 3D shapes other than pure-frequency (pseudo-3D, mixed
# kinetic/freq, etc.) — supported in a later iteration.
function NMRTools.nmrplot(spec::NMRData{<:Any,3}; kwargs...)
    throw(ArgumentError("3D nmrplot in MakieExt currently only supports \
                         pure-frequency 3D spectra. For pseudo-3D, project \
                         to 2D first (e.g. `maximum(spec; dims=2)[:, 1, :]`) \
                         and call nmrplot on the result."))
end

end # module
