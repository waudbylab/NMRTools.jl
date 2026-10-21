# ─────────────────────────────────────────────────────────────────────────
# 1D frequency-domain spectra → lines
# ─────────────────────────────────────────────────────────────────────────

const _Spec1DFreq = NMRData{T,1,<:Tuple{<:FrequencyDimension}} where {T}
const _Spec1DNonFreq = NMRData{T,1,<:Tuple{<:NonFrequencyDimension}} where {T}

# Axis defaults for a 1D frequency spectrum: a single clean bottom axis —
# all left/right/top spines and the entire y-decoration are hidden.
function _axis1d_defaults(dim, title_str)
    return (; xreversed=true,
            xlabel=axislabel(dim),
            ylabel="",
            xgridvisible=false,
            ygridvisible=false,
            yticksvisible=false,
            yticklabelsvisible=false,
            leftspinevisible=false,
            rightspinevisible=false,
            topspinevisible=false,
            bottomspinevisible=true,
            xtickalign=1,
            title=title_str)
end

# Apply x/y limits to a 1D axis. When only `xlims` is given, the y-range is
# auto-scaled to the data *within* the x-window (so an out-of-window peak
# doesn't flatten the view). `offsets` are per-spectrum vstack offsets.
function _apply_1d_limits!(ax, specs, normalize, xlims, ylims; offsets=nothing)
    if !isnothing(xlims) && isnothing(ylims)
        _autoscale_to_xlims!(ax, specs, normalize, xlims; offsets)
    elseif !isnothing(xlims) || !isnothing(ylims)
        _apply_axis_limits!(ax, xlims, ylims)
    end
    return ax
end

function NMRTools.nmrplot(spec::_Spec1DFreq; figure=NamedTuple(), kwargs...)
    fig = Makie.Figure(; figure...)
    ax, plt = NMRTools.nmrplot(fig[1, 1], spec; kwargs...)
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function NMRTools.nmrplot(gp::Union{Makie.GridPosition,Makie.GridSubposition},
                          spec::_Spec1DFreq;
                          normalize=true,
                          xlims=nothing,
                          ylims=nothing,
                          title=nothing,
                          xlabel=nothing,
                          ylabel=nothing,
                          legend=false,
                          axis=NamedTuple(),
                          kwargs...)
    Afwd = reorder(spec, ForwardOrdered)
    dim = dims(Afwd, 1)
    title_str = isempty(refdims(Afwd)) ? string(something(NMRTools.label(Afwd), "")) :
                refdims_title(Afwd)
    ax_kwargs = _axis_overrides(_axis1d_defaults(dim, title_str);
                                title, xlabel, ylabel, axis)
    ax = Makie.Axis(gp; ax_kwargs...)
    plt = NMRTools.nmrplot!(ax, spec; normalize, kwargs...)
    _apply_legend!(ax, legend)
    _apply_1d_limits!(ax, [spec], normalize, xlims, ylims)
    return Makie.AxisPlot(ax, plt)
end

function NMRTools.nmrplot!(ax::Makie.Axis, spec::_Spec1DFreq;
                           normalize=true, kwargs...)
    ax.xreversed = true
    Afwd = reorder(spec, ForwardOrdered)
    sf = _normalization_divisor(Afwd, normalize)
    x = data(dims(Afwd, 1))
    y = _realdata(Afwd) ./ sf
    return Makie.lines!(ax, x, y;
                        label=string(something(NMRTools.label(Afwd), "")),
                        kwargs...)
end

# ─────────────────────────────────────────────────────────────────────────
# 1D non-frequency-domain data → lines for FIDs, scatter otherwise
# ─────────────────────────────────────────────────────────────────────────

function NMRTools.nmrplot(spec::_Spec1DNonFreq; figure=NamedTuple(), kwargs...)
    fig = Makie.Figure(; figure...)
    ax, plt = NMRTools.nmrplot(fig[1, 1], spec; kwargs...)
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function NMRTools.nmrplot(gp::Union{Makie.GridPosition,Makie.GridSubposition},
                          spec::_Spec1DNonFreq;
                          normalize=true,
                          title=nothing,
                          xlabel=nothing,
                          ylabel=nothing,
                          legend=false,
                          axis=NamedTuple(),
                          kwargs...)
    Afwd = reorder(spec, ForwardOrdered)
    dim = dims(Afwd, 1)
    title_str = isempty(refdims(Afwd)) ? string(something(NMRTools.label(Afwd), "")) :
                refdims_title(Afwd)
    defaults = (; xlabel=axislabel(dim),
                ylabel="Intensity",
                xgridvisible=false,
                ygridvisible=false,
                xtickalign=1,
                ytickalign=1,
                title=title_str)
    ax_kwargs = _axis_overrides(defaults; title, xlabel, ylabel, axis)
    ax = Makie.Axis(gp; ax_kwargs...)
    plt = NMRTools.nmrplot!(ax, spec; normalize, kwargs...)
    _apply_legend!(ax, legend)
    return Makie.AxisPlot(ax, plt)
end

function NMRTools.nmrplot!(ax::Makie.Axis, spec::_Spec1DNonFreq;
                           normalize=true, kwargs...)
    Afwd = reorder(spec, ForwardOrdered)
    sf = _normalization_divisor(Afwd, normalize)
    x = data(dims(Afwd, 1))
    y = _realdata(Afwd) ./ sf
    plotfn! = _isfid(dims(Afwd, 1)) ? Makie.lines! : Makie.scatter!
    return plotfn!(ax, x, y;
                   label=string(something(NMRTools.label(Afwd), "")),
                   kwargs...)
end

# a time dimension with an imaginary unit holds an FID, which is drawn as a line
_isfid(dim) = dim isa TimeDimension && !isnothing(get(metadata(dim), :mcindex, nothing))

# ─────────────────────────────────────────────────────────────────────────
# Vector of 1D spectra (overlay or vertical stack)
# ─────────────────────────────────────────────────────────────────────────

function NMRTools.nmrplot(v::AbstractVector{<:NMRData{<:Any,1,
                                                      <:Tuple{<:FrequencyDimension}}};
                          figure=NamedTuple(), kwargs...)
    fig = Makie.Figure(; figure...)
    ax, plt = NMRTools.nmrplot(fig[1, 1], v; kwargs...)
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function NMRTools.nmrplot(gp::Union{Makie.GridPosition,Makie.GridSubposition},
                          v::AbstractVector{<:NMRData{<:Any,1,
                                                      <:Tuple{<:FrequencyDimension}}};
                          normalize=true,
                          vstack=false,
                          xlims=nothing,
                          ylims=nothing,
                          color=nothing,
                          colors=nothing,
                          colormap=nothing,
                          title=nothing,
                          xlabel=nothing,
                          ylabel=nothing,
                          legend=false,
                          axis=NamedTuple(),
                          kwargs...)
    isempty(v) && throw(ArgumentError("nmrplot: empty spectrum vector"))
    dim = dims(reorder(first(v), ForwardOrdered), 1)
    ax_kwargs = _axis_overrides(_axis1d_defaults(dim, ""); title, xlabel, ylabel, axis)
    ax = Makie.Axis(gp; ax_kwargs...)
    colors = isnothing(colors) ? color : colors
    plt = _plot_1d_series!(ax, v; normalize, vstack, colors, colormap, kwargs...)
    _apply_legend!(ax, legend)
    vdelta = _vstack_delta(v, normalize, vstack)
    offsets = [(i - 1) * vdelta for i in eachindex(v)]
    _apply_1d_limits!(ax, v, normalize, xlims, ylims; offsets)
    return Makie.AxisPlot(ax, plt)
end

function _plot_1d_series!(ax::Makie.Axis, v;
                          normalize=true, vstack=false,
                          colors=nothing, colormap=nothing, kwargs...)
    n = length(v)
    cs = _series_poscolors(colors, n; colormap)
    vdelta = _vstack_delta(v, normalize, vstack)

    local first_plt
    voffset = 0.0
    for (i, A) in enumerate(v)
        Afwd = reorder(A, ForwardOrdered)
        sf = _normalization_divisor(Afwd, normalize)
        x = data(dims(Afwd, 1))
        y = _realdata(Afwd) ./ sf .+ voffset
        plt = Makie.lines!(ax, x, y; color=cs[i],
                           label=string(something(NMRTools.label(Afwd), "")),
                           kwargs...)
        i == 1 && (first_plt = plt)
        voffset += vdelta
    end
    return first_plt
end
