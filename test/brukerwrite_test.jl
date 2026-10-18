using NMRTools
using MulticomplexNumbers
using Test

const TESTDATA = joinpath(@__DIR__, "test-data")

"raw on-disk values of a fid/ser file, `td` points per row, padding removed"
function rawvalues(dir, nrows)
    aq = loadjdx(joinpath(dir, "acqus"))
    f = isfile(joinpath(dir, "ser")) ? joinpath(dir, "ser") : joinpath(dir, "fid")
    T = get(aq, :dtypa, 0) == 2 ? Float64 : Int32
    x = Vector{T}(undef, filesize(f) ÷ sizeof(T))
    read!(f, x)
    x = Float64.(get(aq, :bytorda, 0) == 1 ? ntoh.(x) : ltoh.(x))
    return reshape(x, :, nrows)[1:aq[:td], :]
end

@testset "savefid: round trip $name" for (name, nrows) in (("15n-trosy-ts3", 192),
                                                         ("19f-cest-ts3", 101))
    src = joinpath(TESTDATA, name)
    if !isfile(joinpath(src, "ser"))
        @test_skip "raw test data not available"
    else
        fid = loadfid(src)
        mktempdir() do dir
            out = joinpath(dir, "1")
            file = savefid(out, fid)
            @test basename(file) == "ser"
            @test_throws ArgumentError savefid(out, fid)          # no silent overwrite
            savefid(out, fid; overwrite=true)

            # loadfid recovers the data, and the raw values match the original file
            fid2 = loadfid(out)
            @test size(fid2) == size(fid)
            @test parent(fid2) ≈ parent(fid)
            @test rawvalues(out, nrows) ≈ rawvalues(src, nrows)
            # parameters keep their types (TopSpin needs integer TD)
            @test loadjdx(joinpath(out, "acqus"))[:td] isa Integer
            @test loadjdx(joinpath(out, "acqus"))[:dtypa] == 2

            # Int32 output (older TopSpin)
            savefid(out, fid; dtype=Int32, overwrite=true)
            @test loadjdx(joinpath(out, "acqus"))[:dtypa] == 0
            @test rawvalues(out, nrows) ≈ rawvalues(src, nrows)
        end
    end
end

@testset "savefid: processed time-domain data" begin
    src = joinpath(TESTDATA, "15n-trosy-ts3")
    if !isfile(joinpath(src, "ser"))
        @test_skip "raw test data not available"
    else
        fid = loadfid(src) |> linearpredict(T2Dim; npoints=128)
        mktempdir() do dir
            savefid(dir, fid)
            fid2 = loadfid(dir)
            @test size(fid2) == (1536, 128)
            @test parent(fid2) ≈ parent(fid)
        end
        @test_throws ArgumentError savefid(tempname(), ft(fid))
    end
end

@testset "savenmr: round trip" begin
    src = joinpath(TESTDATA, "15n-trosy-ts3")
    if !isfile(joinpath(src, "ser"))
        @test_skip "raw test data not available"
    else
        spec = loadfid(src) |> apodize(Cos2Window()) |> zerofill |> ft |> phase(1, 30.0)
        mktempdir() do dir
            pdir = savenmr(dir, spec)
            @test pdir == joinpath(dir, "pdata", "1")
            @test all(isfile, joinpath.(pdir, ["2rr", "2ir", "2ri", "2ii", "procs", "proc2s"]))
            @test isfile(joinpath(dir, "acqus"))                  # copied from the source
            @test_throws ArgumentError savenmr(dir, spec)

            spec2 = loadnmr(pdir; allcomponents=true)
            @test size(spec2) == size(spec)
            @test collect(data(dims(spec2, 1))) ≈ collect(data(dims(spec, 1)))
            @test collect(data(dims(spec2, 2))) ≈ collect(data(dims(spec, 2)))
            scale = maximum(abs, realest.(parent(spec)))
            @test maximum(abs.(parent(spec2) .- parent(spec))) < 1e-6 * scale
            @test loadjdx(joinpath(pdir, "procs"))[:phc0] ≈ 30.0

            # real part only, into a second procno
            savenmr(dir, real(spec); procno=2)
            @test !isfile(joinpath(dir, "pdata", "2", "2ii"))
            @test data(loadnmr(joinpath(dir, "pdata", "2"))) ≈ data(real(spec)) rtol = 1e-6
        end
        @test_throws ArgumentError savenmr(tempname(), loadfid(src))
    end
end

@testset "savenmr: generated parameter files" begin
    # a spectrum with no source experiment: minimal acqus files are written
    n = 256
    t = collect((0:(n - 1)) ./ 1000.0)
    dd = Dict{Symbol,Any}(:mcindex => 1, :bf => 500.0e6, :sf => 500.0e6, :swhz => 1000.0,
                          :swppm => 2.0, :offsethz => 0.0, :offsetppm => 4.7, :td => n,
                          :label => "1H", :units => "s")
    fid = NMRData([Multicomplex(cos(2π * 100 * ti), sin(2π * 100 * ti)) for ti in t],
                  (T1Dim(t; metadata=dd),))
    spec = ft(fid)
    mktempdir() do dir
        savenmr(dir, spec)
        @test isfile(joinpath(dir, "acqus"))
        spec2 = loadnmr(joinpath(dir, "pdata", "1"); allcomponents=true)
        @test collect(data(dims(spec2, 1))) ≈ collect(data(dims(spec, 1)))
        @test real.(data(spec2)) ≈ realest.(parent(spec)) rtol = 1e-6
    end
    mktempdir() do dir
        savefid(dir, fid)
        fid2 = loadfid(dir)
        # no first-point factor was recorded, so loadfid applies its default scaling of 0.5
        @test parent(fid2)[2:end] ≈ parent(fid)[2:end]
        @test parent(fid2)[1] ≈ 0.5 * parent(fid)[1]
    end
end
