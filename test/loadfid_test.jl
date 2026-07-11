using NMRTools
using MulticomplexNumbers
using Test

const TESTDATA = joinpath(@__DIR__, "test-data")

@testset "loadfid: 2D echo-antiecho (15N TROSY)" begin
    dir = joinpath(TESTDATA, "15n-trosy-ts3")
    if !isfile(joinpath(dir, "ser"))
        @test_skip "raw test data not available"
    else
        dat = loadfid(joinpath(dir, "ser"))

        @test dat isa NMRData
        @test dat[:ndim] == 2
        @test dat[:format] == :bruker_fid

        # phase-sensitive 2D → bicomplex element type
        @test eltype(parent(dat)) <: Multicomplex
        @test order(first(parent(dat))) == 2

        # direct TD=3072 → 1536 complex; indirect TD=192 (EA) → 96 complex
        @test size(dat) == (1536, 96)

        @test dims(dat, 1) isa T1Dim
        @test dims(dat, 2) isa T2Dim
        @test metadata(dat, 1, :mcindex) == 1
        @test metadata(dat, 2, :mcindex) == 2
        @test metadata(dat, 2, :quadrature) == :echoantiecho

        # time axes start at zero and step by the dwell time
        @test data(dims(dat, 1))[1] == 0.0
        @test data(dims(dat, 1))[2] ≈ 1 / metadata(dat, 1, :swhz)

        # routing through loadnmr gives the same result
        dat2 = loadnmr(joinpath(dir, "ser"))
        @test parent(dat2) ≈ parent(dat)
    end
end

@testset "loadfid: pseudo-2D (19F CEST)" begin
    dir = joinpath(TESTDATA, "19f-cest-ts3")
    if !isfile(joinpath(dir, "ser"))
        @test_skip "raw test data not available"
    else
        dat = loadfid(joinpath(dir, "ser"))

        @test dat isa NMRData
        @test dat[:ndim] == 2
        # non-phase-sensitive indirect dimension → order-1 multicomplex
        @test order(first(parent(dat))) == 1
        @test dims(dat, 1) isa T1Dim
        @test metadata(dat, 1, :mcindex) == 1
        # direct complex points are self-consistent with acqus TD
        @test size(dat, 1) == acqus(dat, :td) ÷ 2

        # removefilter=false should also load without error
        raw = loadfid(joinpath(dir, "ser"); removefilter=false)
        @test size(raw) == size(dat)
    end
end

@testset "loadfid: process to spectrum" begin
    dir = joinpath(TESTDATA, "15n-trosy-ts3")
    if !isfile(joinpath(dir, "ser"))
        @test_skip "raw test data not available"
    else
        dat = loadfid(joinpath(dir, "ser"))
        spec = dat |> apodize |> zerofill |> ft |> real
        @test dims(spec, 1) isa FrequencyDimension
        @test dims(spec, 2) isa FrequencyDimension
        @test eltype(parent(spec)) == Float64
        # zero-filled by the default factor of two in each Fourier dimension
        @test size(spec) == (2 * 1536, 2 * 96)
    end
end
