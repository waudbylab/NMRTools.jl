@compile_workload begin
    # 1D data
    axH = F1Dim(8:0.1:9; metadata=Dict{Symbol,Any}(:test => 123))
    dat = NMRData(0.0:1:10, (axH,))

    foo = dat[At(8.5)]
    foo = dat[Near(8.5)]
    foo = dat[8.4 .. 8.6]
    foo = dat[8.41 .. 8.6]
    foo = isnothing(metadata(dat, 1)[:window])
    foo = label(dims(dat, 1))
    foo = isnothing(units(dims(dat, 1)))
    foo = metadata(dat, 1)[:test]

    # 2D data
    x = 8:0.1:9
    y = 100:120
    z = x .+ y'  # 11 x 21

    axH = F1Dim(x)
    axN = F2Dim(y)
    dat = NMRData(z, (axH, axN))

    newax = TkinDim(110:130)
    dat = replacedimension(dat, 2, newax)

    # data import, using the small Bruker datasets bundled with the tests
    # (the example data are lazy artifacts, which would need a download)
    testdata = joinpath(@__DIR__, "..", "test", "test-data")
    if isdir(testdata)
        with_logger(NullLogger()) do
            # 1D, with acqus, title and sample metadata
            spec = loadnmr(joinpath(testdata, "samples", "2"))
            foo = spec[:noise]
            foo = acqus(spec, :pulprog)
            foo = spec[1, :window]
            foo = label(spec)
            foo = sample(spec, :sample, :label)
            foo = spec[8.0 .. 9.0]
            foo = ppm(dims(spec, 1))
            foo = hz(dims(spec, 1))
            foo = detectpeaks(spec)

            # pseudo-2D, with pulse programme annotations
            spec = loadnmr(joinpath(testdata, "19f-cest-ts3"))
            foo = annotations(spec, "cest.offset")
            foo = spec[:, 1]
            foo = detectpeaks(spec[:, 1])

            # raw 1D data: the processing pipeline and automatic processing
            fid = loadfid(joinpath(testdata, "samples", "2"))
            spec = fid |> apodize(ExponentialWindow(1.0)) |> zerofill |> ft |> phase(F1Dim, 10.0)
            foo = realest(spec)
            foo = autoprocess(fid)

            # raw 2D data (bicomplex)
            fid = loadfid(joinpath(testdata, "15n-trosy-ts3"))
            foo = fid |> apodize(Cos2Window()) |> zerofill |> ft |> phase(F2Dim, 0.0, 0.0) |> realest
        end
    end
end
