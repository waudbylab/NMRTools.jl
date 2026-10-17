# ─────────────────────────────────────────────────────────────────────────
# Pole plots: shift × linewidth, area ∝ integral, colour = doubled phase
# ─────────────────────────────────────────────────────────────────────────

const _POLE_COLORSCHEME = ColorSchemes.romaO     # cyclic, so ±180° meet

# colour of a pole: doubled phase mapped onto the cyclic scheme (2φ = 0 at the centre)
function _polecolor(p::Pole)
    t = (rad2deg(angle(cis(2 * angle(p)))) + 180) / 360
    return get(_POLE_COLORSCHEME, t)
end

function NMRTools.poleplot(spec::NMRData; figure=NamedTuple(), kwargs...)
    fig = Makie.Figure(; figure...)
    ax, plt = NMRTools.poleplot(fig[1, 1], spec; kwargs...)
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function NMRTools.poleplot(gp::Union{Makie.GridPosition,Makie.GridSubposition}, spec::NMRData;
                           dim=1,
                           poles=nothing,
                           spectrum=true,
                           colorbar=true,
                           markersize=(4, 24),
                           strokewidth=1.5,
                           title=nothing,
                           axis=NamedTuple(),
                           kwargs...)
    ps = isnothing(poles) ? NMRTools.estimatepoles(spec, dim; kwargs...) : poles
    fdim = dims(spec, dim)
    title_str = something(title, string(something(NMRTools.label(spec), "")))
    layout = Makie.GridLayout(gp)

    showspectrum = spectrum && ndims(spec) == 1
    row = showspectrum ? 2 : 1
    ax = Makie.Axis(layout[row, 1];
                    xreversed=true,
                    yscale=log10,
                    ytickformat=vs -> [string(round(v; sigdigits=2)) for v in vs],
                    xlabel=axislabel(fdim),
                    ylabel="Linewidth (Hz)",
                    xgridvisible=false,
                    ygridvisible=false,
                    xtickalign=1,
                    ytickalign=1,
                    title=showspectrum ? "" : title_str,
                    axis...)
    if showspectrum
        axs = Makie.Axis(layout[1, 1]; _axis1d_defaults(fdim, title_str)...,
                         xlabel="", xticklabelsvisible=false)
        NMRTools.nmrplot!(axs, spec; color=:black, linewidth=1)
        Makie.linkxaxes!(axs, ax)
        Makie.rowsize!(layout, 1, Makie.Relative(0.35))
        Makie.rowgap!(layout, 1, 4)
    end

    plt = NMRTools.poleplot!(ax, ps; markersize, strokewidth)
    if colorbar
        Makie.Colorbar(layout[row, 2];
                       colormap=_POLE_COLORSCHEME.colors,
                       limits=(-180, 180),
                       ticks=-180:90:180,
                       label="2 × phase (°)")
    end
    return Makie.AxisPlot(ax, plt)
end

function NMRTools.poleplot!(ax::Makie.Axis, poles::AbstractVector{Pole};
                            markersize=(4, 24), strokewidth=1.5)
    ax.xreversed = true
    # draw large markers first so that small ones stay visible on top
    ps = sort(poles; by=p -> abs(p.amplitude), rev=true)
    amax = isempty(ps) ? 1.0 : abs(first(ps).amplitude)
    msmin, msmax = markersize
    x = [p.shift for p in ps]
    y = [p.linewidth for p in ps]
    ms = [msmin + (msmax - msmin) * sqrt(abs(p.amplitude) / amax) for p in ps]   # area ∝ integral
    c = [_polecolor(p) for p in ps]
    pos = [cos(angle(p)) >= 0 for p in ps]
    neg = .!pos

    plt = Makie.scatter!(ax, x[pos], y[pos]; color=c[pos], markersize=ms[pos], strokewidth=0)
    Makie.scatter!(ax, x[neg], y[neg]; color=:transparent, strokecolor=c[neg],
                   markersize=ms[neg], strokewidth=strokewidth)
    return plt
end
