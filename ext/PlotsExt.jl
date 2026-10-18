module PlotsExt

"""
notes on testing:
]activate --temp
]dev .
]add Plots
using NMRTools, Plots
"""

using NMRTools
using Plots
using SimpleTraits
using Colors
using ColorSchemes

struct ContourLike end

contourlevels(spacing=1.7, n=12) = (spacing^i for i in 0:(n - 1))

_parse_colorant(c::Colorant) = c
_parse_colorant(c) = parse(Colorant, c)

function _derive_negcolor(poscolor)
    hsv = convert(HSV, _parse_colorant(poscolor))
    HSV(hsv.h, hsv.s * 0.4, min(1.0, max(0.5, hsv.v + 0.4)))
end

axislabel(dat::NMRData, n=1) = axislabel(dims(dat, n))
axislabel(dim::FrequencyDimension) = "$(label(dim)) chemical shift (ppm)"
function axislabel(dim::NMRDimension)
    if isnothing(units(dim))
        "$(label(dim))"
    else
        "$(label(dim)) ($(units(dim)))"
    end
end

_realdata(A::AbstractNMRData) = data(A)
_realdata(A::AbstractNMRData{<:Multicomplex}) = realest.(data(A))

function _checknormalize(normalize)
    normalize isa Bool || normalize isa AbstractNMRData ||
        throw(ArgumentError("normalize must be true, false or a reference spectrum"))
    return nothing
end

"""
    _scaling(A, normalize)

Divisor applied to the intensities of `A` when plotting: `scale(A)` (scans × receiver gain ×
concentration) when normalising, so that spectra recorded under different conditions can be
compared, and 1 otherwise.
"""
function _scaling(A, normalize)
    _checknormalize(normalize)
    return normalize === false ? 1 : scale(A)
end

"""
    _contourbase(A, normalize, ref=A)

Base contour level of `A`, in its raw intensity units; contours are drawn at multiples of this.
Without normalisation it is the noise of `A`. With normalisation it is the noise of the reference
spectrum (`normalize` itself if a spectrum is given, otherwise `ref`) converted to the intensity
scale of `A`, so equal contour levels mean equal normalised intensity across spectra. For a single
spectrum that is its own reference the two coincide.
"""
function _contourbase(A, normalize, ref=A)
    _checknormalize(normalize)
    normalize === false && return A[:noise]
    refspec = normalize isa AbstractNMRData ? normalize : ref
    return refspec[:noise] * scale(A) / scale(refspec)
end

# 1D plot (frequency)
@recipe function f(A::NMRData{T,1,Tuple{D}}; normalize=true) where {T,D<:FrequencyDimension}
    Afwd = reorder(A, ForwardOrdered) # make sure data axes are in forwards order
    x = dims(Afwd, 1)

    # recommend 1D to be a line plot
    seriestype --> :path
    markershape --> :none

    # set default title
    title --> ifelse(isempty(refdims(Afwd)), label(Afwd), refdims_title(Afwd))

    # turn off legend, but provide a label in case the user wishes
    legend --> false
    label --> label(A)

    xguide --> axislabel(A)
    xflip --> true
    xgrid --> false
    xtick_direction --> :out

    yguide --> ""
    yshowaxis --> false
    yticks --> nothing

    scaling = _scaling(Afwd, normalize)

    delete!(plotattributes, :normalize)
    return data(x), _realdata(Afwd) ./ scaling
end

# 1D plot (non-frequency)
@recipe function f(A::NMRData{T,1,Tuple{D}};
                   normalize=true) where {T,D<:NonFrequencyDimension}
    Afwd = reorder(A, ForwardOrdered) # make sure data axes are in forwards order
    x = dims(Afwd, 1)

    # recommend 1D to be a scatter plot
    seriestype --> :scatter

    # set default title
    title --> ifelse(isempty(refdims(Afwd)), label(Afwd), refdims_title(Afwd))

    # turn off legend, but provide a label in case the user wishes
    legend --> false
    label --> label(A)

    xguide --> axislabel(A)
    xflip --> false
    xtick_direction --> :out

    yguide --> "Intensity"
    yshowaxis --> true
    yticks --> :out
    # yerror --> A[:noise] * ones(length(A)) ./ (normalize ? scale(Afwd) : 1)

    widen --> true
    if minimum(x) < 0
        if maximum(x) < 0
            xlims --> (-Inf, 0)
        end
    else
        xlims --> (0, Inf)
    end
    if minimum(A) < 0
        if maximum(A) < 0
            ylims --> (-Inf, 0)
        end
    else
        ylims --> (0, Inf)
    end

    grid --> false
    frame --> :box

    scaling = _scaling(Afwd, normalize)
    delete!(plotattributes, :normalize)
    return data(x), _realdata(Afwd) ./ scaling
end

# multiple 1D plots
@recipe function f(v::Vector{<:NMRData{T,1} where {T}}; normalize=true, vstack=false)
    # recommend 1D to be a line plot
    seriestype --> :path
    markershape --> :none

    # use the first entry to determine axis label
    xguide --> axislabel(v[1])
    xflip --> true
    xgrid --> false
    xtick_direction --> :out

    yguide --> ""
    yshowaxis --> false
    yticks --> nothing

    delete!(plotattributes, :vstack)
    delete!(plotattributes, :normalize)

    voffset = 0
    if vstack isa Bool
        if vstack
            vdelta = maximum([maximum(abs.(A)) / _scaling(A, normalize)
                              for A in v]) / length(v)
            ordering = :back
        else
            vdelta = 0
            ordering = :front
        end
    elseif vstack isa Number
        vdelta = maximum([maximum(abs.(A)) / _scaling(A, normalize)
                          for A in v]) / length(v) * vstack
        ordering = :back
    else
        throw(ArgumentError("vstack must be a Bool or Number"))
    end

    # TODO add guide lines
    # if vstack
    #     yticks --> voffset .+ (0:length(v)-1)*vdelta
    # else
    #     yticks --> [0,]
    # end

    for (i, A) in enumerate(v)
        scaling = _scaling(A, normalize)
        @series begin
            seriestype --> :path
            markershape --> :none
            z_order --> ordering
            Afwd = reorder(A, ForwardOrdered) # make sure data axes are in forwards order
            x = dims(Afwd, 1)
            label --> label(A)
            data(x), _realdata(Afwd) ./ scaling .+ voffset
        end
        voffset += vdelta
    end
end

# 2D plot
@recipe function f(d::D; normalize=true,
                   usegradient=false,
                   poscolor=nothing,
                   negcolor=nothing,
                   negcontours=true) where {D<:NMRData{T,2} where {T}}
    plotattributes[:poscolor] = poscolor
    plotattributes[:negcolor] = negcolor
    plotattributes[:negcontours] = negcontours
    return SimpleTraits.trait(HasNonFrequencyDimension{D}), d
end

@recipe function f(::Type{Not{HasNonFrequencyDimension{D}}},
                   d::D) where {D<:NMRData{T,2} where {T}}
    seriestype --> :contour

    dfwd = reorder(d, ForwardOrdered) # make sure data axes are in forwards order
    # dfwd = DimensionalData.maybe_permute(dfwd, (YDim, XDim))
    x, y = dims(dfwd)

    # set default title
    title --> label(d)
    legend --> false
    framestyle --> :box

    xguide --> axislabel(x)
    xflip --> true
    xgrid --> false
    xtick_direction --> :out

    yguide --> axislabel(y)
    yflip --> true
    ygrid --> false
    ytick_direction --> :out

    normalize = get(plotattributes, :normalize, true)
    delete!(plotattributes, :normalize)

    σ = _contourbase(dfwd, normalize)

    stype = get(plotattributes, :seriestype, nothing)
    if stype ∈ [:heatmap, :wireframe]
        # heatmap
        data(x), data(y), permutedims(_realdata(dfwd))
    else
        poscolor = get(plotattributes, :poscolor, nothing)
        negcolor_arg = get(plotattributes, :negcolor, nothing)
        negcontours = get(plotattributes, :negcontours, true)
        delete!(plotattributes, :poscolor)
        delete!(plotattributes, :negcolor)
        delete!(plotattributes, :negcontours)

        if isnothing(poscolor)
            poscolor = get(plotattributes, :seriescolor, nothing)
        end
        if isnothing(poscolor)
            # Advance through Plots palette for sequential plot!/plot calls.
            # Each spectrum contributes 2 series (pos+neg), so series_count÷2 gives spectrum index.
            p = get(plotattributes, :plot_object, nothing)
            if isnothing(p)
                # Plots.current() throws if no plot exists yet (first plot() call)
                try
                    p = Plots.current()
                catch
                    p = nothing
                end
            end
            n_prev = (p isa Plots.Plot) ? length(p.series_list) : 0
            pal = Plots.palette(:auto)
            poscolor = pal[mod1(n_prev ÷ 2 + 1, length(pal))]
        end

        poscolor_c = _parse_colorant(poscolor)
        negcolor_c = isnothing(negcolor_arg) ? _derive_negcolor(poscolor_c) :
                     _parse_colorant(negcolor_arg)

        @series begin
            levels --> 5σ .* contourlevels()
            seriescolor := poscolor_c
            primary --> true
            data(x), data(y), permutedims(_realdata(dfwd))
        end
        if negcontours
            @series begin
                levels --> -5σ .* contourlevels()
                seriescolor := negcolor_c
                primary := false
                data(x), data(y), permutedims(_realdata(dfwd))
            end
        end
    end
end

# pseudo-2D
# TODO FIX need to treat differently if user asks for a heatmap
@recipe function f(::Type{HasNonFrequencyDimension{D}},
                   d::D) where {D<:NMRData{T,2} where {T}}
    z = reorder(d, ForwardOrdered) # make sure data axes are in forwards order
    # z = DimensionalData.maybe_permute(z, (YDim, XDim))
    if dims(z)[1] isa NonFrequencyDimension
        z = transpose(z)
    end
    x, y = dims(z)

    # set default title
    title --> label(d)
    legend --> false
    framestyle --> :box

    xguide --> axislabel(x)
    xflip --> true
    xgrid --> false
    xtick_direction --> :out

    yguide --> axislabel(y)
    yflip --> false
    ygrid --> false
    ytick_direction --> :out

    normalize = get(plotattributes, :normalize, true)
    delete!(plotattributes, :normalize)

    scaling = _scaling(z, normalize)

    stype = get(plotattributes, :seriestype, nothing)
    if stype ∈ [:heatmap, :wireframe]
        # heatmap
        data(x), data(y), permutedims(_realdata(z)) ./ scaling
    else
        # default

        # TODO - BUG - this attribute isn't recognised
        gradient = get(plotattributes, :usegradient, false)
        delete!(plotattributes, :usegradient)

        palettesize = length(y)
        if palettesize > 64
            palettesize = floor(palettesize / 4)
        elseif palettesize > 48
            palettesize = floor(palettesize / 3)
        elseif palettesize > 32
            palettesize = floor(palettesize / 2)
        elseif palettesize > 16
            palettesize = 16
        end

        # don't set a gradient palette if colours already specified
        setpalette = :seriescolor ∉ keys(plotattributes) ||
                     :linecolor ∉ keys(plotattributes)

        zreal = _realdata(z)
        xones = ones(length(x))
        for i in 1:length(y)
            @series begin
                seriestype --> :path3d
                seriescolor --> i

                if gradient
                    line_z --> zreal[:, i]  # colour by height
                    palette --> :darkrainbow
                else
                    if setpalette
                        palette --> palette(:phase, 11)
                    end
                end
                primary --> (i == 1)
                fillrange --> 0 # not currently implemented in plots for 3D data
                data(x), xones * y[i], zreal[:, i] / scaling
            end
        end
    end
end

# multiple 2D plots
@recipe function f(v::Vector{D}; normalize=true,
                   poscolor=nothing,
                   negcolor=nothing,
                   negcontours=true) where {D<:NMRData{T,2}} where {T}
    plotattributes[:poscolor] = poscolor
    plotattributes[:negcolor] = negcolor
    plotattributes[:negcontours] = negcontours
    return SimpleTraits.trait(HasNonFrequencyDimension{D}), v
end

@recipe function f(::Type{Not{HasNonFrequencyDimension{D}}},
                   v::Vector{D}) where {D<:NMRData{T,2}} where {T}
    n = length(v)

    # Resolve poscolors / negcolors
    poscolors_arg = get(plotattributes, :poscolor, nothing)
    negcolors_arg = get(plotattributes, :negcolor, nothing)
    negcontours = get(plotattributes, :negcontours, true)
    delete!(plotattributes, :poscolor)
    delete!(plotattributes, :negcolor)
    delete!(plotattributes, :negcontours)

    if isnothing(poscolors_arg)
        poscolors_arg = get(plotattributes, :seriescolor, nothing)
    end

    if isnothing(poscolors_arg)
        poscolors = [HSV(h, 0.9, 0.85) for h in (0:(n - 1)) .* (360.0 / n)]
    elseif isa(poscolors_arg, AbstractVector)
        poscolors = [_parse_colorant(poscolors_arg[mod1(i, length(poscolors_arg))]) for i in 1:n]
    else
        poscolors = fill(_parse_colorant(poscolors_arg), n)
    end

    if isnothing(negcolors_arg)
        negcolors = [_derive_negcolor(c) for c in poscolors]
    elseif isa(negcolors_arg, AbstractVector)
        negcolors = [_parse_colorant(negcolors_arg[mod1(i, length(negcolors_arg))]) for i in 1:n]
    else
        negcolors = fill(_parse_colorant(negcolors_arg), n)
    end

    seriestype --> :contour

    dfwd = reorder(v[1], ForwardOrdered) # make sure data axes are in forwards order
    x, y = dims(dfwd)

    # set default title
    title --> ""
    colorbar --> nothing
    framestyle --> :box

    xguide --> axislabel(x)
    xflip --> true
    xgrid --> false
    xtick_direction --> :out

    yguide --> axislabel(y)
    yflip --> true
    ygrid --> false
    ytick_direction --> :out

    normalize = get(plotattributes, :normalize, true)
    delete!(plotattributes, :normalize)

    # contour levels follow each spectrum's own noise (normalize=false), or the noise of the first
    # spectrum or a given reference converted to each spectrum's scale

    for (i, d) in enumerate(v)
        dfwd = reorder(d, ForwardOrdered) # make sure data axes are in forwards order
        x, y = dims(dfwd)
        poscolor_c = poscolors[i]
        negcolor_c = negcolors[i]

        σ = _contourbase(dfwd, normalize, first(v))

        @series begin
            levels --> 5σ .* contourlevels()
            seriescolor := poscolor_c
            primary := false
            label := nothing
            data(x), data(y), permutedims(_realdata(dfwd))
        end
        if negcontours
            @series begin
                levels --> -5σ .* contourlevels()
                seriescolor := negcolor_c
                primary := false
                label := nothing
                data(x), data(y), permutedims(_realdata(dfwd))
            end
        end
        @series begin
            seriestype := :path
            seriescolor := poscolor_c
            primary --> true
            label --> label(dfwd)
            [], []
        end
    end
end

@recipe function f(::Type{HasNonFrequencyDimension{D}},
                   v::Vector{D}) where {D<:NMRData{T,2}} where {T}
    # Axis setup from first dataset
    z0 = reorder(v[1], ForwardOrdered)
    if dims(z0)[1] isa NonFrequencyDimension
        z0 = transpose(z0)
    end
    x0, y0 = dims(z0)

    title --> ""
    legend --> false
    framestyle --> :box

    xguide --> axislabel(x0)
    xflip --> true
    xgrid --> false
    xtick_direction --> :out

    yguide --> axislabel(y0)
    yflip --> false
    ygrid --> false
    ytick_direction --> :out

    normalize = get(plotattributes, :normalize, true)
    delete!(plotattributes, :normalize)

    for (j, d) in enumerate(v)
        z = reorder(d, ForwardOrdered)
        if dims(z)[1] isa NonFrequencyDimension
            z = transpose(z)
        end
        x, y = dims(z)

        scaling = _scaling(z, normalize)

        zreal = _realdata(z)
        xones = ones(length(x))
        for i in 1:length(y)
            @series begin
                seriestype --> :path3d
                seriescolor --> j
                primary --> (i == 1)
                fillrange --> 0
                data(x), xones * y[i], zreal[:, i] / scaling
            end
        end
    end
end

# Plot 1D peaks as vertical lines
@recipe function f(peaks::AbstractVector{<:Peak{1}})
    seriestype --> :vline
    linecolor --> :red
    linestyle --> :solid
    alpha --> 0.3
    label --> nothing
    z_order --> 1

    return [ppm(p, 1) for p in peaks]
end

# Plot 2D peaks as crosses
@recipe function f(peaks::AbstractVector{<:Peak{2}})
    seriestype --> :scatter
    markershape --> :xcross
    markercolor --> :black
    markersize --> 3
    label --> nothing
    xflip --> true
    yflip --> true

    return [ppm(p, 1) for p in peaks], [ppm(p, 2) for p in peaks]
end

end # module