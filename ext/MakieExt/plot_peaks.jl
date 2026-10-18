# ─────────────────────────────────────────────────────────────────────────
# Peak lists → overlays on 1D and 2D spectra
# ─────────────────────────────────────────────────────────────────────────

# 1D peaks are drawn as vertical lines, with their labels along the top of the axis.
# Lines need no y-coordinate, so they sit correctly over a spectrum however it was
# normalised or stacked, and the labels follow the axis as it zooms.
function NMRTools.nmrplot!(ax::Makie.Axis, pks::AbstractVector{<:Peak{1}};
                           color=:red, alpha=0.4, labels=true, labelsize=10,
                           kwargs...)
    x = [ppm(p, 1) for p in pks]
    plt = Makie.vlines!(ax, x; color=(color, alpha), kwargs...)
    if labels && any(p -> !isempty(NMRTools.label(p)), pks)
        pos = Makie.lift(ax.finallimits) do lims
            ytop = lims.origin[2] + lims.widths[2]
            return [Makie.Point2f(xi, ytop) for xi in x]
        end
        Makie.text!(ax, pos; text=NMRTools.label.(pks), color, fontsize=labelsize,
                    align=(:center, :top), offset=(0, -2))
    end
    return plt
end

# 2D peaks are drawn as crosses, labelled just above and to the right.
function NMRTools.nmrplot!(ax::Makie.Axis, pks::AbstractVector{<:Peak{2}};
                           color=:black, marker=:xcross, markersize=8, labels=true,
                           labelsize=10, kwargs...)
    ax.xreversed = true
    ax.yreversed = true
    pos = [Makie.Point2f(ppm(p)...) for p in pks]
    plt = Makie.scatter!(ax, pos; color, marker, markersize, kwargs...)
    if labels && any(p -> !isempty(NMRTools.label(p)), pks)
        Makie.text!(ax, pos; text=NMRTools.label.(pks), color, fontsize=labelsize,
                    align=(:left, :bottom), offset=(3, 3))
    end
    return plt
end

function NMRTools.nmrplot!(pks::AbstractVector{<:Peak}; kwargs...)
    ax = Makie.current_axis()
    isnothing(ax) && throw(ArgumentError("nmrplot!: no current axis. \
        Call nmrplot(spec) first, or pass an Axis explicitly."))
    return NMRTools.nmrplot!(ax, pks; kwargs...)
end

# A peak list on its own, without a spectrum: e.g. a list predicted from BMRB shifts. The
# axes are set up as for a spectrum, but peaks carry no nucleus, so the axis labels are
# generic unless given.
function NMRTools.nmrplot(pks::AbstractVector{<:Peak}; figure=NamedTuple(), kwargs...)
    fig = Makie.Figure(; figure...)
    ap = NMRTools.nmrplot(fig[1, 1], pks; kwargs...)
    return Makie.FigureAxisPlot(fig, ap.axis, ap.plot)
end

function NMRTools.nmrplot(gp::Union{Makie.GridPosition,Makie.GridSubposition},
                          pks::AbstractVector{<:Peak{1}};
                          title="", xlabel="Chemical shift (ppm)", axis=NamedTuple(),
                          kwargs...)
    defaults = (; xreversed=true, xlabel, title, xgridvisible=false, ygridvisible=false,
                yticksvisible=false, yticklabelsvisible=false, leftspinevisible=false,
                rightspinevisible=false, topspinevisible=false, xtickalign=1)
    ax = Makie.Axis(gp; merge(defaults, axis)...)
    plt = NMRTools.nmrplot!(ax, pks; kwargs...)
    Makie.ylims!(ax, 0, 1)
    return Makie.AxisPlot(ax, plt)
end

function NMRTools.nmrplot(gp::Union{Makie.GridPosition,Makie.GridSubposition},
                          pks::AbstractVector{<:Peak{2}};
                          title="", xlabel="Chemical shift (ppm)",
                          ylabel="Chemical shift (ppm)", axis=NamedTuple(), kwargs...)
    defaults = (; xreversed=true, yreversed=true, xlabel, ylabel, title,
                xgridvisible=false, ygridvisible=false, xtickalign=1, ytickalign=1)
    ax = Makie.Axis(gp; merge(defaults, axis)...)
    plt = NMRTools.nmrplot!(ax, pks; kwargs...)
    return Makie.AxisPlot(ax, plt)
end
