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
    @test apodize(spec, T1Dim, Cos²Window()) isa NMRData
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

@testset "processing: real/imag/abs/complex" begin
    spec = synthetic_fid()
    r = real(spec)
    @test eltype(parent(r)) == Float64
    @test parent(r) ≈ realest.(parent(spec))

    @test eltype(parent(abs(spec))) == Float64
    @test eltype(parent(complex(spec))) <: Complex
    # complex conversion preserves the components
    @test real.(parent(complex(spec))) ≈ realest.(parent(spec))
end

@testset "processing: pipeline and collections" begin
    spec = synthetic_fid()
    out = spec |> apodize |> zerofill |> ft |> real
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
end
