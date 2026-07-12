using NMRTools
using MulticomplexNumbers
using Test

const NB = NMRTools.NMRBase

# Build a synthetic 1D time-domain NMRData (order-1 multicomplex) from a complex FID.
function fid_from_complex(z::AbstractVector{<:Complex}; sw=1000.0, bf=500.0e6)
    n = length(z)
    t = collect((0:(n - 1)) ./ sw)
    dd = Dict{Symbol,Any}(:mcindex => 1, :bf => bf, :sf => bf, :swhz => sw,
                          :swppm => 1e6 * sw / bf, :offsethz => 0.0, :offsetppm => 0.0,
                          :td => n, :tdzf => n, :npoints => n, :region => missing,
                          :window => nothing, :pseudodim => false, :units => "s")
    d = T1Dim(t; metadata=dd)
    return NMRData([Multicomplex(real(zi), imag(zi)) for zi in z], (d,))
end

@testset "harmonic core: pole recovery" begin
    N = 64
    ω, λ = 2π * 0.11, 0.02          # per-sample angular frequency and damping
    x = ComplexF64[exp((im * ω - λ) * n) for n in 0:(N - 1)]
    poles = NB._hankelpoles(x, 1)
    @test poles[1] ≈ exp(im * ω - λ) rtol = 1e-8

    # two components, both recovered
    ω2, λ2 = -2π * 0.2, 0.05
    x2 = ComplexF64[exp((im * ω - λ) * n) + exp((im * ω2 - λ2) * n) for n in 0:(N - 1)]
    p2 = NB._hankelpoles(x2, 2)
    for zt in (exp(im * ω - λ), exp(im * ω2 - λ2))
        @test minimum(abs.(p2 .- zt)) < 1e-6
    end
end

@testset "harmonic core: linear prediction (coefficient method)" begin
    N = 48
    z1 = exp(im * 2π * 0.07 - 0.03)                       # one damped exponential
    x = ComplexF64[z1^n for n in 0:(N - 1)]
    ext = NB._lppredict(x; npred=16, order=1, direction=:forward)
    @test length(ext) == N + 16
    @test ext ≈ ComplexF64[z1^n for n in 0:(N + 16 - 1)] rtol = 1e-6

    z2 = exp(-im * 2π * 0.2 - 0.05)                       # two components → order 2
    x2 = ComplexF64[z1^n + z2^n for n in 0:(N - 1)]
    ext2 = NB._lppredict(x2; npred=12, order=2, direction=:forward)
    @test ext2 ≈ ComplexF64[z1^n + z2^n for n in 0:(N + 12 - 1)] rtol = 1e-5

    # first-point repair reconstructs early points by backward prediction
    rep = NB._lprepair(x; nrepair=4, order=1)
    @test length(rep) == N
    @test rep ≈ x rtol = 1e-5
end

@testset "processing: frequencyshift" begin
    z = ComplexF64[exp(2π * im * 100 * (n / 1000)) for n in 0:63]   # 100 Hz signal
    spec = fid_from_complex(z; sw=1000.0)
    # shift there and back is the identity
    back = frequencyshift(frequencyshift(spec, T1Dim, 100.0), T1Dim, -100.0)
    @test parent(back) ≈ parent(spec)
    # after moving 100 Hz onto the carrier, the FID is (near) constant → small variation
    shifted = frequencyshift(spec, T1Dim, 100.0)
    v = realest.(parent(shifted))
    @test maximum(v) - minimum(v) < 1e-6
end

@testset "processing: linearpredict (multicomplex glue)" begin
    N = 40
    ω, λ = 2π * 0.09, 0.02
    z = ComplexF64[exp((im * ω - λ) * n) for n in 0:(N - 1)]
    spec = fid_from_complex(z)
    lp = linearpredict(spec, T1Dim; factor=2, order=1)
    @test size(lp, 1) == 2N
    @test metadata(lp, 1, :td) == 2N
    # predicted region matches the analytic continuation
    got = [realest(m) + imag(m) * im for m in parent(lp)]   # back to Complex (order-1)
    analytic = ComplexF64[exp((im * ω - λ) * n) for n in 0:(2N - 1)]
    @test got ≈ analytic rtol = 1e-5
end

@testset "processing: solventfilter lowpass" begin
    n, sw = 128, 1000.0
    t = (0:(n - 1)) ./ sw
    solvent = 50.0 .* exp.(-t .* 1.0)                       # on-resonance (0 Hz), strong, slow
    peak = exp.(2π * im * 250 .* t) .* exp.(-t .* 8.0)      # off-resonance
    spec = fid_from_complex(ComplexF64.(solvent) .+ peak; sw=sw)

    filt = solventfilter(spec; width=40.0, shape=:sine)
    @test size(filt) == size(spec)
    # the strong on-resonance (low-frequency) component is suppressed
    before = abs(sum(realest.(parent(spec))))
    after = abs(sum(realest.(parent(filt))))
    @test after < 0.2 * before
end
