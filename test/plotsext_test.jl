using NMRTools
using Artifacts
using LazyArtifacts
using Plots
using Test
using VisualRegressionTests

"""
notes on testing, to generate plots:

]activate --temp
]dev .
]add Plots
using NMRTools, Plots
"""

@testset "PlotsExt: 1D 19F" begin
    dat = exampledata("1D_19F")

    pr = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(), dat)

    @test pr[1].plotattributes[:xguide] == "19F chemical shift (ppm)"
    @test pr[1].plotattributes[:title] == "Example 19F 1D spectrum"
    @test pr[1].plotattributes[:xflip] == true
    @test pr[1].plotattributes[:xtick_direction] == :out
    @test pr[1].plotattributes[:yshowaxis] == false

    @test pr[1].args[1][1] ≈ -124.16252020246354
    @test pr[1].args[2][1] ≈ 1.2039054352746015
end

@testset "PlotsExt: 2D HN" begin
    dat = exampledata("2D_HN")

    pr = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(), dat)

    @test pr[1].plotattributes[:normalize] == true
    @test pr[1].args[2][1, 1] ≈ 20.85546875

    @test label(pr[1].args[2]) == "13C,15N ubiquitin"
end

# @testset "PlotsExt: normalize parameter with reference spectrum" begin
#     dat1 = exampledata("2D_HN")
#     dat2 = exampledata("2D_HN") * 2  # scaled version

#     # Test single 2D spectrum with reference normalization
#     pr_ref = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => dat1), dat2)
#     pr_self = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => true), dat2)
#     pr_none = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => false), dat2)

#     # When using reference normalization, the scaling should be based on dat1
#     @test pr_ref[1].plotattributes[:levels] ≠ pr_self[1].plotattributes[:levels]
#     @test pr_ref[1].plotattributes[:levels] ≠ pr_none[1].plotattributes[:levels]

#     # Test vector of 2D spectra with reference normalization
#     dat_vector = [dat1, dat2]
#     pr_vector_ref = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => dat1), dat_vector)
#     pr_vector_true = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => true), dat_vector)

#     # Both should use consistent contour levels when using reference spectrum
#     @test length(pr_vector_ref) >= 6  # Should have multiple series for contours and legends
#     @test length(pr_vector_true) >= 6

#     # Test 1D spectrum with reference normalization
#     dat1_1d = exampledata("1D_19F")
#     dat2_1d = dat1_1d * 3  # scaled version

#     pr_1d_ref = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => dat1_1d), dat2_1d)
#     pr_1d_self = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => true), dat2_1d)

#     # Data should be scaled differently when using reference vs self normalization
#     @test pr_1d_ref[1].args[2] ≠ pr_1d_self[1].args[2]

#     # Test error handling
#     @test_throws Exception Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(:normalize => "invalid"), dat1)
# end

@testset "PlotsExt: Multicomplex 2D" begin
    datMC = loadnmr(artifact"2D_HN"; allcomponents=true)
    @test parent(datMC) isa Matrix{<:Multicomplex}

    pr = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(), datMC)

    @test pr[1].plotattributes[:normalize] == true
    @test label(pr[1].args[2]) == "13C,15N ubiquitin"
end

@testset "PlotsExt: Multicomplex 1D (slice)" begin
    datMC2D = loadnmr(artifact"2D_HN"; allcomponents=true)
    datMC = datMC2D[:, 1]
    @test parent(datMC) isa AbstractVector{<:Multicomplex}

    # 1D recipe goes all the way to series data in one apply_recipe call
    pr = Plots.RecipesBase.apply_recipe(Dict{Symbol,Any}(), datMC)
    @test pr[1].plotattributes[:xflip] == true
    @test eltype(pr[1].args[2]) <: Real
end

@testset "PlotsExt: visual regression testing" begin
    dat1 = exampledata("1D_19F")
    @plottest begin
        plot(dat1)
    end "plots-images/plot_1D_19F.png" false 0.02

    dat2 = exampledata("2D_HN")
    @plottest begin
        plot(dat2)
    end "plots-images/plot_2D_HN.png" false 0.02
    @plottest begin
        plot(dat2; c=:black)
    end "plots-images/plot_2D_black.png" false 0.02
    @plottest begin
        plot(dat2; poscolor=:red)
    end "plots-images/plot_2D_red.png" false 0.02
    @plottest begin
        plot(dat2; poscolor=:blue, negcolor=:red)
    end "plots-images/plot_2D_blue_red.png" false 0.02
    @plottest begin
        plot(dat2; poscolor=:limegreen, negcontours=false)
    end "plots-images/plot_2D_noneg.png" false 0.02

    dat3 = exampledata("pseudo2D_XSTE")
    dat3 = setgradientlist(dat3, LinRange(0.05, 0.95, 10), 0.55)
    @plottest begin
        plot(dat3[7 .. 9, :])
    end "plots-images/plot_pseudo2D_XSTE.png" false 0.02
    @plottest begin
        heatmap(dat3[7 .. 9, :])
    end "plots-images/plot_pseudo2D_XSTE_heatmap.png" false 0.02

    dat4 = exampledata("1D_19F_titration")
    @plottest begin
        plot(dat4; xlims=(-125, -122))
    end "plots-images/plot_1D_19F_titration.png" false 0.02
    @plottest begin
        plot(dat4; xlims=(-125, -121), vstack=true)
    end "plots-images/plot_1D_19F_titration_vstack.png" false 0.02
    @plottest begin
        plot(dat4; xlims=(-125, -121), vstack=5)
    end "plots-images/plot_1D_19F_titration_vstack_5.png" false 0.02

    dat5 = exampledata("2D_HN_titration")
    @plottest begin
        plot(dat5)
    end "plots-images/plot_2D_HN_titration.png" false 0.02
    @plottest begin
        plot(dat5[1])
        plot!(dat5[11]; c=:hotpink)
        xlims!(7, 8)
        ylims!(115, 120)
    end "plots-images/plot_2D_HN_overlay.png" false 0.02
    @plottest begin
        plot(dat5[1])
        plot!(dat5[11]; normalize=dat5[1])
        xlims!(7, 8)
        ylims!(115, 120)
    end "plots-images/plot_2D_overlay_normalised.png" false 0.02
    @plottest begin
        plot(dat5[1]; poscolor=:black, negcolor=:red)
        plot!(dat5[11]; normalize=dat5[1], poscolor=:blue, negcolor=:orange)
        xlims!(6.5, 7.5)
        ylims!(114, 118)
    end "plots-images/plot_2D_overlay_colors.png" false 0.02

    @plottest begin
        dat6 = exampledata("3D_HNCA")
        plot(maximum(dat6; dims=2)[:, 1, :] / 3)
    end "plots-images/plot_3D_HNCA_projection.png" false 0.02
end
