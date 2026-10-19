using NMRTools
using MulticomplexNumbers
using Test

# Build a synthetic 1D time-domain NMRData (order-1 multicomplex FID) with a single frequency.
function synthetic_fid(; n=32, sw=1000.0, bf=500.0e6, ν=100.0)
    t = collect((0:(n - 1)) ./ sw)
    dd = Dict{Symbol,Any}(:mcindex => 1,
                          :bf => bf,
                          :sf => bf,
                          :swhz => sw,
                          :swppm => 1e6 * sw / bf,
                          :offsethz => 0.0,
                          :offsetppm => 0.0,
                          :td => n,
                          :tdzf => n,
                          :npoints => n,
                          :region => missing,
                          :window => nothing,
                          :pseudodim => false,
                          :units => "s")
    d = T1Dim(t; metadata=dd)
    y = [Multicomplex(cos(2π * ν * ti), sin(2π * ν * ti)) for ti in t]
    return NMRData(y, (d,))
end

@testset "processing: apodize" begin
    spec = synthetic_fid()
    a = apodize(spec, T1Dim, ExponentialWindow(10.0))
    @test a isa NMRData
    @test size(a) == size(spec)
    @test eltype(parent(a)) == eltype(parent(spec))
    # exponential apodization attenuates later points
    @test abs(parent(a)[end]) < abs(parent(spec)[end])
    # window is recorded in metadata
    @test metadata(a, 1, :window) isa ExponentialWindow

    # inverse apodization undoes it
    b = apodize(a, T1Dim, ExponentialWindow(10.0); inverse=true)
    @test parent(b) ≈ parent(spec)

    # all-dimensions form with an explicit window
    @test size(apodize(spec, ExponentialWindow(10.0))) == size(spec)

    # no-argument form apodizes the (only) Fourier dimension
    @test size(apodize(spec)) == size(spec)

    # a window with no length constructs and gets tmax filled in automatically
    @test apodize(spec, T1Dim, Cos2Window()) isa NMRData
end

@testset "processing: zerofill" begin
    spec = synthetic_fid(; n=32)
    z = zerofill(spec, T1Dim; factor=2)
    @test size(z, 1) == 64
    @test metadata(z, 1, :tdzf) == 64
    # original points preserved, tail zeroed
    @test parent(z)[1:32] ≈ parent(spec)
    @test all(iszero, parent(z)[33:64])
    # default factor doubles
    @test size(zerofill(spec)) == (64,)
end

@testset "processing: ft round-trip" begin
    spec = synthetic_fid()
    f = ft(spec, T1Dim)
    @test dims(f, 1) isa FrequencyDimension
    @test size(f) == size(spec)
    @test units(f, 1) == "ppm"

    # inverse FT recovers the FID
    back = ft(f, F1Dim; inverse=true)
    @test dims(back, 1) isa TimeDimension
    @test parent(back) ≈ parent(spec)
end

@testset "processing: phase" begin
    spec = synthetic_fid()
    f = ft(spec, T1Dim)
    # a full 360° zero-order phase is the identity
    @test parent(phase(f, F1Dim, 360.0)) ≈ parent(f)
    # zero phase is the identity
    @test parent(phase(f, F1Dim, 0.0)) ≈ parent(f)
    # curried form
    p = phase(F1Dim, 90.0)
    @test parent(p(f)) ≈ parent(phase(f, F1Dim, 90.0))
end

@testset "processing: realest/real/imag/abs/complex" begin
    spec = synthetic_fid()
    r = realest(spec)
    @test r isa NMRData
    @test eltype(parent(r)) == Float64
    @test parent(r) ≈ realest.(parent(spec))
    # real and imag act elementwise, as in MulticomplexNumbers; for 1D real agrees with realest
    @test parent(real(spec)) ≈ parent(r)
    @test parent(imag(spec)) ≈ imag.(parent(spec))
    @test imag(spec) isa NMRData

    @test eltype(parent(abs(spec))) == Float64
    @test eltype(parent(complex(spec))) <: Complex
    # complex conversion preserves the components
    @test real.(parent(complex(spec))) ≈ realest.(parent(spec))
end

@testset "processing: pipeline and collections" begin
    spec = synthetic_fid()
    out = spec |> apodize |> zerofill |> ft |> realest
    @test eltype(parent(out)) == Float64
    @test dims(out, 1) isa FrequencyDimension
    @test size(out) == (64,)

    # dictionary mapping
    d = Dict(:a => synthetic_fid(), :b => synthetic_fid())
    dout = apodize(d)
    @test dout isa AbstractDict
    @test Set(keys(dout)) == Set([:a, :b])
    @test all(v isa NMRData for v in values(dout))
end

@testset "processing: extract" begin
    spec = synthetic_fid()
    f = ft(spec, T1Dim)
    ax = data(dims(f, 1))
    lo, hi = minimum(ax), maximum(ax)
    mid = (lo + hi) / 2
    sub = extract(F1Dim, mid .. hi)(f)
    @test sub isa NMRData
    @test size(sub, 1) < size(f, 1)
    @test parent(sub) == parent(f[F1Dim(mid .. hi)])
    @test parent(f |> extract(F1Dim => mid .. hi)) == parent(sub)
end

@testset "processing: digital-filter group delay in ft" begin
    # an FID whose t = 0 is g points into the array (filter build-up before it) transforms to the
    # same spectrum as the FID starting at t = 0, followed by g points of zeros
    n, g, sw = 256, 12, 1000.0
    base = synthetic_fid(; n=n, sw=sw)
    x = parent(base)
    delayed = vcat(zeros(eltype(x), g), x[1:(n - g)])
    t = collect(((0:(n - 1)) .- g) ./ sw)
    md = Dict{Symbol,Any}(pairs(metadata(dims(base, 1)).val))
    md[:grpdly] = Float64(g)
    fid = NMRData(delayed, (T1Dim(t; metadata=md),))
    ref = NMRData(vcat(x[1:(n - g)], zeros(eltype(x), g)), (dims(base, 1),))
    @test parent(ft(fid)) ≈ parent(ft(ref))
    @test metadata(ft(fid), 1, :grpdly) == 0
end

@testset "processing: realest/real/imag of 2D data" begin
    testdata = joinpath(@__DIR__, "test-data", "15n-trosy-ts3")
    if !isfile(joinpath(testdata, "ser"))
        @test_skip "raw test data not available"
    else
        spec = loadfid(testdata) |> ft
        m = parent(spec)
        # components of a bicomplex number are (rr, ir, ri, ii); real and imag strip im2
        @test parent(realest(spec)) == component.(m, 1)
        @test parent(real(spec)) == real.(m)
        @test parent(imag(spec)) == imag.(m)
        @test eltype(parent(real(spec))) <: Multicomplex{Float64,1}
        @test realest(real(spec)) == realest(spec)
    end
end

@testset "processing: pipeline forms of analysis functions" begin
    fid = synthetic_fid(; n=512, ν=100.0)
    fid = NMRTools.rebuild(fid; data=parent(fid) .* exp.(-π * 5.0 .* data(dims(fid, 1))))  # 5 Hz line
    clean = ft(zerofill(fid))
    noisy = NMRTools.rebuild(clean; data=parent(clean) .+ [Multicomplex(1e-2 * randn(), 1e-2 * randn())
                                                         for _ in 1:length(clean)])
    spec = estimatenoise(noisy)

    pks = detectpeaks(spec; snr=20)
    @test !isempty(pks)
    @test (spec |> detectpeaks(snr=20)) == pks
    δ = ppm(first(pks), 1)
    @test (spec |> findpeak(δ)) == findpeak(spec, δ)
    @test (pks |> register(spec)) == register(pks, spec)
    @test (copy(pks) |> register!(spec)) == register!(copy(pks), spec)

    @test (spec |> estimatepoles) == estimatepoles(spec)
    @test (spec |> estimatepoles(1; snr=10)) == estimatepoles(spec, 1; snr=10)
    @test (spec |> estimatephase(order=0)).p0 == estimatephase(spec; order=0).p0

    @test data(dims(spec |> shiftdim(F1Dim, 0.1), 1)) == data(dims(shiftdim(spec, F1Dim, 0.1), 1))
    refd = spec |> reference(F1Dim, δ => δ + 0.05; indirect=false)
    @test data(dims(refd, 1)) == data(dims(reference(spec, F1Dim, δ => δ + 0.05; indirect=false), 1))
    @test parent(spec |> decimate(2)) == parent(decimate(spec, 2))

    pseudo = NMRData(randn(8, 3), (F1Dim(range(8.0, 7.0; length=8)), X2Dim(1:3)))
    t = [0.1, 0.2, 0.3]
    @test dims(pseudo |> setrelaxtimes(t, "s"), 2) isa TrelaxDim
    @test data(dims(pseudo |> setrelaxtimes(2, t), 2)) == data(dims(setrelaxtimes(pseudo, 2, t), 2))
    @test dims(pseudo |> setkinetictimes(t), 2) isa TkinDim
    @test dims(pseudo |> replacedimension(2, TrelaxDim(t)), 2) isa TrelaxDim

    mktempdir() do dir
        file = joinpath(dir, "peaklist.csv")
        pks |> savepeaks(file)
        @test length(loadpeaks(file)) == length(pks)
        spec |> savenmr(joinpath(dir, "1"))
        @test isfile(joinpath(dir, "1", "pdata", "1", "1r"))
    end
end

