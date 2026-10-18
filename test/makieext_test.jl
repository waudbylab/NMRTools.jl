using NMRTools
using Artifacts
using LazyArtifacts
using Test
using CairoMakie
using FileIO
using VisualRegressionTests

"""
notes on testing, to generate plots:

]activate --temp
]dev .
]add CairoMakie
using NMRTools, CairoMakie

ENV["VISUAL_REGRESSION_TESTS_AUTO"] = true;
pkg> test Plots
"""

"""
    save_rgb(filename, fig)

Extracts the image data from a Makie Figure, strips the alpha channel 
by converting it to pure RGB, and saves it to disk. 
"""
function save_rgb(filename, fig)
    rgba_buf = CairoMakie.colorbuffer(fig)

    # Use CairoMakie's internal Colors module to access RGB safely
    rgb_buf = CairoMakie.Colors.RGB.(rgba_buf)

    FileIO.save(filename, rgb_buf)
    return nothing
end

@testset "MakieExt: 1D 19F" begin
    dat = exampledata("1D_19F")
    fig, ax, plt = nmrplot(dat)

    @test fig isa Makie.Figure
    @test ax isa Makie.Axis
    @test plt isa Makie.Lines
    @test ax.xreversed[] == true
    @test ax.xlabel[] == "19F chemical shift (ppm)"
    @test ax.xgridvisible[] == false
    @test ax.yticksvisible[] == false

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F.png"
end

@testset "MakieExt: 2D HN" begin
    dat = exampledata("2D_HN")
    fig, ax, plt = nmrplot(dat)

    @test fig isa Makie.Figure
    @test ax isa Makie.Axis
    @test ax.xreversed[] == true
    @test ax.yreversed[] == true
    @test ax.title[] == "13C,15N ubiquitin"
    # Positive + negative contour layers added.
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN.png"
end

@testset "MakieExt: 2D HN colour overrides" begin
    dat = exampledata("2D_HN")
    fig, ax, plt = nmrplot(dat; poscolor=:red)
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_poscolor.png"

    fig, ax, plt = nmrplot(dat; poscolor=:blue, negcolor=:red)
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_negcolor.png"

    fig, ax, plt = nmrplot(dat; poscolor=:limegreen, negcontours=false)
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 1

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_negcontours.png"

    # `color=` is an alias for `poscolor=` — negative contour is auto-derived,
    # not coloured the same as positive.
    fig, ax, plt = nmrplot(dat; color=:red)
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_color.png"

    # Explicit negcolor wins even when color/poscolor is set.
    fig, ax, plt = nmrplot(dat; color=:green, negcolor=:red)
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_explicit_negcolor.png"
end

@testset "MakieExt: top-level axis overrides" begin
    dat = exampledata("2D_HN")
    fig, ax, plt = nmrplot(dat; title="My HN", xlabel="¹H", ylabel="¹⁵N")
    @test ax.title[] == "My HN"
    @test ax.xlabel[] == "¹H"
    @test ax.ylabel[] == "¹⁵N"

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_overrides.png"
end

@testset "MakieExt: legend kwarg" begin
    dats = exampledata("2D_HN_titration")
    fig, ax, plt = nmrplot(dats; legend=true)
    @test any(c -> c isa Makie.Legend, fig.content)

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_legend.png"

    fig, ax, plt = nmrplot(dats; legend=:lt)
    @test any(c -> c isa Makie.Legend, fig.content)

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_legend_topleft.png"

    fig, ax, plt = nmrplot(dats)
    @test !any(c -> c isa Makie.Legend, fig.content)

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_legend_none.png"
end

@testset "MakieExt: nmrplot! forwards from Figure / FigureAxisPlot" begin
    dat = exampledata("2D_HN")
    fap = nmrplot(dat)
    # Both forms should work.
    plt2 = nmrplot!(fap, dat; poscolor=:red)
    @test plt2 isa Makie.Contour
    plt3 = nmrplot!(fap.figure, dat; poscolor=:blue)
    @test plt3 isa Makie.Contour

    @visualtest (fname -> save_rgb(fname, fap)) "makie-images/2D_HN_nmrplotbang.png"
end

@testset "MakieExt: 2D overlay cycles colours" begin
    dat = exampledata("2D_HN")
    fap = nmrplot(dat)
    plt1_color = fap.plot.color[]
    plt2 = nmrplot!(fap, dat / 2)
    plt2_color = plt2.color[]
    # Second spectrum should get a different Wong colour, not the same one.
    @test plt1_color != plt2_color

    @visualtest (fname -> save_rgb(fname, fap)) "makie-images/2D_HN_overlay.png"
end

@testset "MakieExt: nmrplot! into existing axis" begin
    dat = exampledata("1D_19F")
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    plt = nmrplot!(ax, dat)
    @test ax.xreversed[] == true
    @test plt isa Makie.Lines

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_nmrplotbang.png"
end

@testset "MakieExt: grid position" begin
    dat = exampledata("1D_19F")
    fig = Makie.Figure()
    ax, plt = nmrplot(fig[1, 1], dat)
    @test ax isa Makie.Axis
    @test ax.xreversed[] == true

    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_gridpos.png"

    fig = Figure(; size=(1000, 420))
    nmrplot(fig[1, 1], exampledata("1D_19F"); title="1D")
    nmrplot(fig[1, 2], exampledata("2D_HN"); title="2D", xlims=(6, 10))
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/subplots.png"
end

@testset "MakieExt: vector of 1D" begin
    dats = exampledata("1D_19F_titration")
    fig, ax, plt = nmrplot(dats)
    @test ax isa Makie.Axis
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == length(dats)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_titration.png"

    fig, ax, plt = nmrplot(dats; vstack=true)
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == length(dats)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_titration_vstack.png"

    fig, ax, plt = nmrplot(dats; vstack=5)
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == length(dats)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_titration_vstack5.png"
end

@testset "MakieExt: 1D xlims auto-scales y to window" begin
    dat = exampledata("1D_19F")
    # No xlims: limits stay automatic.
    fig, ax, plt = nmrplot(dat)
    @test ax.limits[] === nothing || all(isnothing, ax.limits[])
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_nolims.png"

    # With xlims: explicit limits applied (x window + windowed y range).
    fig, ax, plt = nmrplot(dat; xlims=(-123, -122))
    @test ax.limits[] !== nothing
    fl = ax.finallimits[]
    @test Float64(minimum(fl)[1]) ≈ -123.0 atol = 1e-6
    @test Float64(maximum(fl)[1]) ≈ -122.0 atol = 1e-6
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_xlim.png"

    # Works for a series too.
    dats = exampledata("1D_19F_titration")
    fig, ax, plt = nmrplot(dats; xlims=(-125, -121))
    @test ax.limits[] !== nothing
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_19F_titration_xlim.png"
end

@testset "MakieExt: vector of 2D" begin
    dats = exampledata("2D_HN_titration")
    fig, ax, plt = nmrplot(dats)
    @test ax isa Makie.Axis
    @test ax.xreversed[] == true
    @test ax.yreversed[] == true
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2 * length(dats)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_titration.png"

    # Explicit colours
    fig, ax, plt = nmrplot(dats; colors=[:red, :blue])
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2 * length(dats)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_titration_colors.png"

    # Wong colours for small number of points
    fig, ax, plt = nmrplot(dats[1:5])
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_titration_wong.png"

    # Sequential colormap fallback for long series.
    fig, ax, plt = nmrplot(dats; colormap=:plasma)
    @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2 * length(dats)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_titration_colormap.png"
end

@testset "MakieExt: Multicomplex 1D slice" begin
    datMC2D = loadnmr(artifact"2D_HN"; allcomponents=true)
    datMC = datMC2D[:, 1]
    @test parent(datMC) isa AbstractVector{<:Multicomplex}

    fig, ax, plt = nmrplot(datMC)
    @test ax.xreversed[] == true
    @test plt isa Makie.Lines
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/1D_MC_slice.png"
end

@testset "MakieExt: 2D projections" begin
    dat = exampledata("2D_HN")
    fig, ax, plt = nmrplot(dat; xprojection=true, yprojection=true)
    @test ax isa Makie.Axis
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_projections.png"

    fig, ax, plt = nmrplot(dat; xprojection=:sum, yprojection=:sum)
    @test ax isa Makie.Axis
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_projections_sum.png"

    fig, ax, plt = nmrplot(dat; xprojection=true)
    @test ax isa Makie.Axis
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_xprojection.png"

    slice = dat[:, 1]
    fig, ax, plt = nmrplot(dat; xprojection=slice)
    @test ax isa Makie.Axis
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_xprojection_slice.png"
end

@testset "MakieExt: overlay with projections" begin
    dat = exampledata("2D_HN")
    fap = nmrplot(dat; xprojection=true, yprojection=true)
    plt2 = nmrplot!(fap, dat / 2; xprojection=true, yprojection=true,
                    poscolor=:red)
    @test plt2 isa Makie.Contour
    @visualtest (fname -> save_rgb(fname, fap)) "makie-images/2D_HN_overlay_projections.png"
end

@testset "MakieExt: vector of pseudo-2D" begin
    dat = exampledata("pseudo2D_XSTE")
    fig, ax, plt = nmrplot([dat, dat / 2]; style=:flat)
    @test ax isa Makie.Axis
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_flat_series.png"

    fig, ax, plt = nmrplot([dat, dat / 2]; style=:waterfall)
    @test ax isa Makie.Axis3
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_waterfall_series.png"

    @test_throws ArgumentError nmrplot([dat, dat / 2]; style=:heatmap)
end

@testset "MakieExt: pseudo-2D overlay" begin
    dat = exampledata("pseudo2D_XSTE")
    fap = nmrplot(dat; style=:flat)
    plt2 = nmrplot!(dat / 2; style=:flat, color=:green)
    @test plt2 isa Makie.Lines
    @visualtest (fname -> save_rgb(fname, fap)) "makie-images/pseudo2D_flat_overlay.png"

    fap = nmrplot(dat; style=:waterfall)
    plt3 = nmrplot!(dat / 2; style=:waterfall, color=:green)
    @test plt3 isa Makie.Lines
    # increase tolerence - different makie versions move the axis text very slightly
    @visualtest (fname -> save_rgb(fname, fap)) "makie-images/pseudo2D_waterfall_overlay.png" false 0.1
end

@testset "MakieExt: pseudo-2D heatmap (default)" begin
    dat = exampledata("pseudo2D_XSTE")
    fig, ax, plt = nmrplot(dat)
    @test ax isa Makie.Axis
    @test ax.xreversed[] == true
    @test plt isa Makie.Heatmap
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_heatmap.png"

    # No colorbar
    fig, ax, plt = nmrplot(dat; colorbar=false)
    @test plt isa Makie.Heatmap
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_heatmap_nocolorbar.png"
end

@testset "MakieExt: pseudo-2D flat" begin
    dat = exampledata("pseudo2D_XSTE")
    fig, ax, plt = nmrplot(dat; style=:flat)
    @test ax isa Makie.Axis
    @test ax.xreversed[] == true
    @test plt isa Makie.Lines
    # One Lines per slice; zero baseline adds an HLines (not a Lines).
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == length(dims(dat, 2))
    @test any(p -> p isa Makie.HLines, ax.scene.plots)
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_flat.png"
end

@testset "MakieExt: pseudo-2D waterfall" begin
    dat = exampledata("pseudo2D_XSTE")
    fig, ax, plt = nmrplot(dat; style=:waterfall)
    @test ax isa Makie.Axis3
    @test plt isa Makie.Lines
    @test count(p -> p isa Makie.Lines, ax.scene.plots) == length(dims(dat, 2))
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_waterfall.png"
end

# @testset "MakieExt: 3D contour" begin
#     dat = exampledata("3D_HNCA")
#     fig, ax, plt = nmrplot(dat)
#     @test ax isa Makie.Axis3
#     @test plt isa Makie.Contour
#     @visualtest (fname->save_rgb(fname, fig)) "makie-images/3D_HNCA.png"
#     # Positive + negative iso-surface.
#     @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2

#     fig, ax, plt = nmrplot(dat; negcontours=false)
#     @test count(p -> p isa Makie.Contour, ax.scene.plots) == 1

#     fig, ax, plt = nmrplot(dat; poscolor=:blue, level=30)
#     @test count(p -> p isa Makie.Contour, ax.scene.plots) == 2
# end

# @testset "MakieExt: 3D contour overlay" begin
#     dat = exampledata("3D_HNCA")
#     # Vector overlay, like 2D contour series.
#     fig, ax, plt = nmrplot([dat, dat / 2])
#     @test ax isa Makie.Axis3
#     @test count(p -> p isa Makie.Contour, ax.scene.plots) == 4  # 2 specs × ±

#     # Manual overlay advances the colour cycle.
#     fap = nmrplot(dat)
#     plt2 = nmrplot!(dat / 2)
#     @test plt2 isa Makie.Contour
# end

@testset "MakieExt: general xlims/ylims" begin
    # 2D
    dat2 = exampledata("2D_HN")
    fig, ax, plt = nmrplot(dat2; xlims=(6, 10), ylims=(110, 125))
    @test ax.limits[] !== nothing
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/2D_HN_lims.png"

    # pseudo-2D heatmap
    datp = exampledata("pseudo2D_XSTE")
    fig, ax, plt = nmrplot(datp; xlims=(5, 10))
    @test ax.limits[] !== nothing
    @visualtest (fname -> save_rgb(fname, fig)) "makie-images/pseudo2D_XSTE_xlim.png"
end
