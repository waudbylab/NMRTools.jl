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

export nmrplot, nmrplot!

include("precompile.jl")

end
