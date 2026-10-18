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

        # the direct time axis steps by the dwell time and starts GRPDLY dwells before t = 0
        sw = metadata(dat, 1, :swhz)
        @test data(dims(dat, 1))[1] ≈ -76 / sw
        @test data(dims(dat, 1))[2] - data(dims(dat, 1))[1] ≈ 1 / sw
        @test data(dims(dat, 2))[1] == 0.0

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
        spec = dat |> apodize |> zerofill |> ft |> realest
        @test dims(spec, 1) isa FrequencyDimension
        @test dims(spec, 2) isa FrequencyDimension
        @test eltype(parent(spec)) == Float64
        # zero-filled by the default factor of two in each Fourier dimension
        @test size(spec) == (2 * 1536, 2 * 96)
    end
end

@testset "loadfid: digital filter and first-point rules" begin
    trosy = joinpath(TESTDATA, "15n-trosy-ts3")     # baseopt (DIGMOD = 3), DSPFVS = 21
    cest = joinpath(TESTDATA, "19f-cest-ts3")       # digital (DIGMOD = 1), DSPFVS = 20
    if !isfile(joinpath(trosy, "ser")) || !isfile(joinpath(cest, "ser"))
        @test_skip "raw test data not available"
    else
        # baseopt: group delay recorded, build-up kept at negative time, t = 0 point not scaled
        dat = loadfid(trosy)
        raw = loadfid(trosy; removefilter=false, firstpoint=1.0)
        @test metadata(dat, 1, :grpdly) == 76
        @test metadata(dat, 1, :firstpoint) == 1.0
        @test metadata(dat, 2, :firstpoint) == 0.5
        @test parent(dat)[:, 2:end] ≈ parent(raw)[:, 2:end]        # direct data unchanged
        @test metadata(raw, 1, :grpdly) == 0
        @test data(dims(raw, 1))[1] == 0.0
        # an explicit first-point factor is applied at t = 0 (point 77)
        dat5 = loadfid(trosy; firstpoint=0.5)
        @test metadata(dat5, 1, :firstpoint) == 0.5
        @test parent(dat5)[77, 2] ≈ 0.5 * parent(raw)[77, 2]       # column 2: indirect t > 0
        @test parent(dat5)[76, 2] ≈ parent(raw)[76, 2]
        # ft removes the delay; the spectrum's dimension has none left
        spec = ft(dat, 1)
        @test metadata(spec, 1, :grpdly) == 0
        @test metadata(ft(spec, 1; inverse=true), 1, :grpdly) == 0
        # the window starts at t = 0: build-up points are not apodized
        ap = apodize(dat, 1, ExponentialWindow(50.0))
        @test parent(ap)[1:76, :] == parent(dat)[1:76, :]
        @test !(parent(ap)[200, 1] ≈ parent(dat)[200, 1])

        # digital (not baseopt): fractional group delay, first point scaled by 0.5
        dat = loadfid(cest)
        @test metadata(dat, 1, :grpdly) ≈ 67.9868469238281
        @test metadata(dat, 1, :firstpoint) == 0.5

        # DSPFVS outside 10–23 (e.g. after cryoproc1d): no group-delay removal, with a warning
        mktempdir() do dir
            cp(cest, joinpath(dir, "1"))
            NMRTools.NMRIO._patchjdxfile(joinpath(dir, "1", "acqus"), ["DSPFVS" => 0])
            mod = @test_logs (:warn, r"DSPFVS = 0") match_mode = :any loadfid(joinpath(dir, "1"))
            @test metadata(mod, 1, :grpdly) == 0
            @test parent(mod) ≈ parent(loadfid(cest; removefilter=false))
        end
    end
end

@testset "loadfid: modified FID time origin (cryoproc1d)" begin
    # cryoproc1d removes the filter build-up plus some points after t = 0; acqu keeps the original TD
    mktempdir() do dir
        write(joinpath(dir, "acqu"), "##TITLE= Parameter file, TopSpin 4.0.0\n##\$TD= 65536\n##END=\n")
        acqus = Dict{Symbol,Any}(:td => 65340, :grpdly => 76.0, :dspfvs => 0, :digmod => 3)
        @test NMRTools.NMRIO._modifiedfidoffset(acqus, dir) == 22
        @test (@test_logs (:warn, r"22 points after") NMRTools.NMRIO._groupdelay(acqus, dir)) == -22
        # without acqu the offset is unknown: no delay
        @test (@test_logs (:warn, r"DSPFVS = 0") NMRTools.NMRIO._groupdelay(acqus, nothing)) == 0
    end
end

@testset "loadfid: intensity scaling by 2^NC" begin
    cest = joinpath(TESTDATA, "19f-cest-ts3")             # NC = -6, DTYPA = 0
    if !isfile(joinpath(cest, "ser"))
        @test_skip "raw test data not available"
    else
        raw = Vector{Int32}(undef, 2)
        open(f -> read!(f, raw), joinpath(cest, "ser"))
        dat = loadfid(cest; removefilter=false, firstpoint=1.0)
        @test realest(parent(dat)[1, 1]) ≈ ltoh(raw[1]) * 2.0^-6
        @test imag(parent(dat)[1, 1]) ≈ ltoh(raw[2]) * 2.0^-6
    end
end
