module NMRTools

using PrecompileTools
using Reexport: Reexport
using MulticomplexNumbers

include("NMRBase/NMRBase.jl")
include("NMRIO/NMRIO.jl")

Reexport.@reexport using .NMRBase
Reexport.@reexport using .NMRIO
Reexport.@reexport using MulticomplexNumbers

"""
    nmrplot(spec; kwargs...)
    nmrplot(gridpos, spec; kwargs...)
    nmrplot!(ax, spec; kwargs...)

Plot an NMR spectrum with conventions appropriate for NMR data (reversed
chemical-shift axes, hidden y-decorations on 1D, ±contours for 2D).

Requires `Makie` (or a Makie backend such as `GLMakie`, `CairoMakie`,
`WGLMakie`) to be loaded.
"""
function nmrplot end

"""
    nmrplot!(ax, spec; kwargs...)

Mutating variant of [`nmrplot`](@ref) — plots into an existing `Axis`,
overriding axis orientation/labels to keep NMR conventions intact.
"""
function nmrplot! end

"""
    poleplot(spec; poles=nothing, spectrum=true, colorbar=true, kwargs...)
    poleplot(gridpos, spec; kwargs...)
    poleplot!(ax, poles; markersize=(4, 24), strokewidth=1.5)

Plot the pole model of a frequency-domain spectrum: one marker per [`Pole`](@ref) at its chemical
shift (x) and linewidth (y, log scale). Marker area is proportional to the line integral. Colour shows
twice the line phase on a cyclic colour map, so positive and negative absorption lines have the same
colour and a correctly phased spectrum is a single colour. A zero-order error shifts every colour; a
first-order error appears as a colour gradient along the shift axis. Filled markers have a positive
real part, open markers a negative one.

`poles` defaults to `estimatepoles(spec; kwargs...)`; other keywords are passed to
[`estimatepoles`](@ref). With `spectrum=true` (1D data only) the real part of the spectrum is drawn
above the pole map on a linked axis.

Requires `Makie` (or a Makie backend) to be loaded.
"""
function poleplot end

"""
    poleplot!(ax, poles; markersize=(4, 24), strokewidth=1.5)

Draw a vector of [`Pole`](@ref)s into an existing axis. See [`poleplot`](@ref).
"""
function poleplot! end

export nmrplot, nmrplot!
export poleplot, poleplot!

include("precompile.jl")

end
