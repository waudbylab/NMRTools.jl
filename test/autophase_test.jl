using NMRTools
using MulticomplexNumbers
using Test

# Synthetic 1D spectrum: a sum of Lorentzian lines (shift in Hz from the carrier, linewidth in Hz,
# complex amplitude), Fourier transformed with the standard pipeline.
function synthetic_spectrum(lines; n=4096, sw=5000.0, bf=500.0e6, noise=0.0)
    t = collect((0:(n - 1)) ./ sw)
    z = zeros(ComplexF64, n)
    for (ν, lw, a) in lines
        z .+= a .* cis.(2π * ν .* t) .* exp.(-π * lw .* t)
    end
    z .+= noise .* (randn(n) .+ im .* randn(n))
    z[1] /= 2                                         # t = 0 point carries half weight
    dd = Dict{Symbol,Any}(:mcindex => 1, :bf => bf, :sf => bf, :swhz => sw,
                          :swppm => 1e6 * sw / bf, :offsethz => 0.0, :offsetppm => 4.7,
                          :td => n, :tdzf => n, :npoints => n, :region => missing,
                          :window => nothing, :pseudodim => false, :units => "s")
    fid = NMRData([Multicomplex(real(zi), imag(zi)) for zi in z], (T1Dim(t; metadata=dd),))
    return fid |> zerofill |> ft
end

const LINES = [(-2000.0, 4.0, 1.0), (-1300.0, 6.0, 0.8), (-600.0, 3.0, 1.2), (-100.0, 5.0, 0.6),
               (400.0, 4.0, 1.0), (900.0, 8.0, 0.9), (1500.0, 3.0, 0.7), (2100.0, 5.0, 1.1)]

@testset "estimatepoles: synthetic lines" begin
    spec = synthetic_spectrum(LINES; noise=1e-3)
    poles = estimatepoles(spec)
    @test issorted([p.shift for p in poles]; rev=true)
    bf = 500.0
    for (ν, lw, a) in LINES
        δ = 4.7 + ν / bf
        p = poles[argmin([abs(q.shift - δ) for q in poles])]
        @test abs(p.shift - δ) < 1e-3
        @test p.linewidth ≈ lw rtol = 0.1
        @test abs(angle(p)) < deg2rad(3)            # absorptive lines have zero phase
        @test p.snr > 100
    end
    # excluded regions are dropped
    poles2 = estimatepoles(spec; exclude=[(4.7 + 380 / bf) .. (4.7 + 420 / bf)])
    @test length(poles2) < length(poles)
    @test_throws ArgumentError estimatepoles(ft(synthetic_spectrum(LINES); inverse=true))
end

@testset "estimatephase / autophase: synthetic" begin
    spec = synthetic_spectrum(LINES; noise=1e-3)
    for (p0, p1) in ((0.0, 0.0), (37.0, 0.0), (-120.0, 25.0), (150.0, -60.0))
        bad = phase(spec, 1, p0, p1)
        est = estimatephase(bad)
        @test abs(rad2deg(angle(cis(deg2rad(est.p0 + p0))))) < 3
        @test abs(est.p1 + p1) < 3
        @test est.residual < 5
        @test est.nlines >= length(LINES)
        @test all(r -> r.relsnr < 0.01, est.rejected)
        @test all(e -> e.used, est.edges)                # clean baseline: edges agree with lines

        est0 = estimatephase(phase(spec, 1, p0, 0.0); order=0)
        @test est0.p1 == 0
        @test abs(rad2deg(angle(cis(deg2rad(est0.p0 + p0))))) < 3
    end

    # a minority of negative lines does not flip the result
    neglines = [(ν, lw, i in (2, 5) ? -a : a) for (i, (ν, lw, a)) in enumerate(LINES)]
    spec2 = synthetic_spectrum(neglines; noise=1e-3)
    est = estimatephase(phase(spec2, 1, 60.0, 0.0))
    @test abs(rad2deg(angle(cis(deg2rad(est.p0 + 60.0))))) < 3

    # edges cannot override the lines: with a baseline offset (first point not halved) they are
    # inconsistent with the line fit and left out
    offset = synthetic_spectrum(LINES; noise=1e-3) |> s -> phase(s, 1, 40.0)
    biased = NMRTools.rebuild(offset; data=parent(offset) .+ 0.05 * maximum(abs, parent(offset)))
    est = estimatephase(biased)
    @test abs(rad2deg(angle(cis(deg2rad(est.p0 + 40.0))))) < 3
    @test est.nlines >= length(LINES)
    # without edge observations
    @test isempty(estimatephase(spec; edges=0).edges)

    # a strong line with a different phase is rejected, reported, and warned about
    oddlines = [(ν, lw, i == 3 ? a * cis(deg2rad(70)) : a) for (i, (ν, lw, a)) in enumerate(LINES)]
    spec3 = synthetic_spectrum(oddlines; noise=1e-3)
    est = estimatephase(spec3; order=0)
    @test abs(rad2deg(angle(cis(deg2rad(est.p0))))) < 3
    @test length(est.rejected) == 1
    @test est.rejected[1].shift ≈ 4.7 - 600 / 500 atol = 1e-3
    @test abs(est.rejected[1].residual) > 60
    @test_logs (:info,) (:warn, r"strong lines left out") autophase(spec3; order=0)

    # pivot: p0 reported at the pivot reproduces the same correction
    bad = phase(spec, 1, -50.0, 40.0)
    e1 = estimatephase(bad)
    e2 = estimatephase(bad; pivot=4.7)
    @test real(phase(bad, 1, e2.p0, e2.p1; pivot=4.7)) ≈ real(phase(bad, 1, e1.p0, e1.p1))

    # autophase applies the correction and records the total phase in the metadata
    fixed = autophase(phase(spec, 1, 30.0, 20.0))
    ph = metadata(fixed, 1, :phase)
    @test abs(ph.p0) < 3 && abs(ph.p1) < 3
    @test minimum(data(real(fixed))) > -0.05 * maximum(data(real(fixed)))
    # curried form
    @test real(phase(spec, 1, 30.0) |> autophase) ≈ real(autophase(phase(spec, 1, 30.0)))

    @test_throws NMRToolsError estimatephase(spec; exclude=[-100 .. 100])
    @test_throws ArgumentError estimatephase(spec; order=2)
end

@testset "autophase0: synthetic" begin
    spec = synthetic_spectrum(LINES; noise=1e-3)
    for p0 in (0.0, 37.0, -120.0, 150.0)
        fixed = autophase0(phase(spec, 1, p0))
        @test abs(rad2deg(angle(cis(deg2rad(metadata(fixed, 1, :phase).p0))))) < 1
        @test metadata(fixed, 1, :phase).p1 == 0
    end
    # excluded regions are left out of the integral
    bad = phase(spec, 1, 40.0)
    @test real(autophase0(bad; exclude=[4.7 + 380 / 500 .. 4.7 + 420 / 500])) isa NMRData
    # curried form
    @test real(bad |> autophase0) ≈ real(autophase0(bad))
    @test_throws ArgumentError autophase0(ft(spec; inverse=true))
end

@testset "phase: metadata and sign convention" begin
    spec = synthetic_spectrum(LINES)
    @test isnothing(metadata(spec, 1, :phase))
    s1 = phase(spec, 1, 30.0, 10.0)
    @test metadata(s1, 1, :phase).p0 ≈ 30.0
    @test metadata(s1, 1, :phase).p1 ≈ 10.0
    @test isnothing(metadata(spec, 1, :phase))           # input unchanged
    s2 = phase(s1, 1, -10.0, 5.0)
    @test metadata(s2, 1, :phase).p0 ≈ 20.0
    @test metadata(s2, 1, :phase).p1 ≈ 15.0
    # a pivot is converted to the first-point convention
    s3 = phase(spec, 1, 30.0, 10.0; pivot=4.7)
    n = length(spec)
    ipiv = NMRTools.NMRBase._pivotindex(dims(spec, 1), 4.7)
    @test metadata(s3, 1, :phase).p0 ≈ 30.0 - 10.0 * ipiv / n
    # integer arguments (previously ambiguous with the curried form)
    @test phase(spec, 1, 0) isa NMRData
    # the TopSpin sign convention: data are multiplied by exp(-im1·θ)
    @test parent(phase(spec, 1, 90.0))[1] ≈ parent(spec)[1] * (-im1)
end
