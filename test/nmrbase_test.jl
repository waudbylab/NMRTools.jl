using NMRTools
using Artifacts
using LazyArtifacts
using Test

@testset "NMRBase: NMRData" begin
    dat1a = exampledata("pseudo2D_XSTE")
    dat1b = exampledata("pseudo2D_XSTE")
    dat2 = exampledata("1D_19F")
    dat2b = NMRData(dat2) # create NMRData from another NMRData object
    @test dat1a == dat1b # test equality
    @test dat1a ≠ dat2
    @test dat2b == dat2

    @test Array(dat1a) isa Matrix{Float64}
    @test Array(dat1a)[1, 1] == 167.84442138671875

    @test collect(dat1a) isa Matrix{Float64}
    @test collect(dat1a)[1, 1] == 167.84442138671875

    dat3 = dat1a[:, 3]  # slicing
    @test ndims(dat3) == 1
    @test dat3 isa NMRData{Float64,1}
    @test size(dat3) == (2048,)

    dat3b = dat1a[100, :]  # slicing
    @test ndims(dat3b) == 1
    @test size(dat3b) == (10,)

    @test data(dat2)[1] == 4671.5184326171875
    @test data(dat2, 1)[1] == -119.5922

    @test hasnonfrequencydimension(dat1a) == true
    @test hasnonfrequencydimension(dat2) == false
end

@testset "NMRBase: dimensions" begin
    axH = F1Dim(8:0.1:9)
    dat = NMRData(0.0:1:10, (axH,))

    @test dims(dat, 1) isa FrequencyDimension
    @test dat[1] == 0
    @test dat[1:3] == [0.0, 1.0, 2.0]
    @test dat[end] == 10
    @test length(dat) == 11
    @test size(dat) == (11,)
    @test (2dat .+ 1)[[1, end]] == [1.0, 21.0]

    @test dat[At(8.5)] == 5
    @test dat[Near(8.5)] == 5
    @test dat[Near(8.49)] == 5
    @test dat[Near(8.54)] == 5
    @test dat[8.4 .. 8.6] == [4, 5, 6]
    @test dat[8.41 .. 8.6] == [5, 6]

    @test_throws DimensionMismatch NMRData(0.0:1:9, (axH,))  # size of axis and data don't match

    # 2D data
    x = 8:0.1:9
    y = 100:120
    z = x .+ y'  # 11 x 21
    axH = F1Dim(x)
    axN = F2Dim(y)
    dat = NMRData(z, (axH, axN))

    @test size(dat) == (11, 21)
    @test length(dat) == 11 * 21
    @test dat[1, :] == y .+ 8
    @test dat[:, 1] == x .+ 100
    @test dat[8.2 .. 8.4, 101 .. 102.5] == [109.2 110.2; 109.3 110.3; 109.4 110.4]

    # pseudo-2D data
    axP = X2Dim(y)
    dat = NMRData(z, (axH, axP))
    datK = setkinetictimes(dat, LinRange(0, 2, 21), "min")
    @test units(datK, 2) == "min"
    datR = setrelaxtimes(dat, LinRange(0, 2, 21))
    @test units(datR, 2) == ""
    @test_logs (:warn,
                "a maximum gradient strength of 0.55 T m⁻¹ is being assumed - this is roughly correct for modern Bruker systems but calibration is recommended") setgradientlist(dat,
                                                                                                                                                                                 LinRange(0.05,
                                                                                                                                                                                          0.95,
                                                                                                                                                                                          21))
    datG = setgradientlist(dat, LinRange(0.05, 0.95, 21), 0.55)
    @test data(datG, 2)[end] == 0.5225
    @test units(datG, 2) == "T m⁻¹"

    # offsets
    spec = exampledata("1D_1H")
    spec2 = shiftdim(spec, 1, 0.5)
    @test data(spec2, 1)[1] - data(spec, 1)[1] == 0.5
    spec2 = shiftdim(spec2, F1Dim, 0.5)
    @test data(spec2, 1)[1] - data(spec, 1)[1] == 1

    # Test ppm() and hz() functions for FrequencyDimension
    ax = dims(spec, 1)

    # Test ppm(axis) - should return same as data(axis)
    @test ppm(ax) == data(ax)
    @test ppm(ax)[1] isa Float64

    # Test hz(axis) - should convert ppm to Hz offsets
    hz_offsets = hz(ax)
    bf = ax[:bf]
    offsetppm = ax[:offsetppm]
    @test hz_offsets[1] ≈ bf * (data(ax)[1] - offsetppm) * 1e-6
    @test hz_offsets[end] ≈ bf * (data(ax)[end] - offsetppm) * 1e-6

    # Test hz(δ, axis) - should convert single ppm value to Hz
    test_ppm = 7.5
    hz_single = hz(test_ppm, ax)
    @test hz_single ≈ bf * (test_ppm - offsetppm) * 1e-6

    # Test hz(δ, axis) with array of values
    test_ppms = [7.0, 7.5, 8.0]
    hz_array = hz(test_ppms, ax)
    @test length(hz_array) == 3
    @test hz_array[1] ≈ bf * (7.0 - offsetppm) * 1e-6
    @test hz_array[2] ≈ bf * (7.5 - offsetppm) * 1e-6
    @test hz_array[3] ≈ bf * (8.0 - offsetppm) * 1e-6
end

@testset "NMRBase: metadata" begin
    # create some test data
    axH = F1Dim(8:0.1:9)
    dat = NMRData(0.0:1:10, (axH,))

    # NMRData metadata
    @test metadata(dat) isa Metadata{NMRData}
    @test metadata(dat) == defaultmetadata(NMRData)
    @test metadata(dat)[:title] == ""
    label!(dat, "my spectrum")
    @test label(dat) == "my spectrum"
    @test isnothing(units(dat))

    # dimension metadata
    axH = F1Dim(8:0.1:9)
    dat = NMRData(0.0:1:10, (axH,))
    @test metadata(dat, 1) isa Metadata{FrequencyDimension}
    @test metadata(dat, 1) == defaultmetadata(FrequencyDimension)
    @test isnothing(metadata(dat, 1)[:window])
    @test label(dims(dat, 1)) == ""
    @test isnothing(units(dims(dat, 1)))

    axH = F1Dim(8:0.1:9; metadata=NoMetadata())
    dat = NMRData(0.0:1:10, (axH,))
    @test metadata(dat, 1) isa Metadata{FrequencyDimension}
    @test metadata(dat, 1) == defaultmetadata(FrequencyDimension)
    @test isnothing(metadata(dat, 1)[:window])
    label!(dat, 1, "my label")
    @test label(dims(dat, 1)) == "my label"
    @test isnothing(units(dims(dat, 1)))

    axH = F1Dim(8:0.1:9; metadata=Dict{Symbol,Any}(:test => 123))
    dat = NMRData(0.0:1:10, (axH,))
    @test metadata(dat, 1) isa Metadata{FrequencyDimension}
    @test issubset(defaultmetadata(FrequencyDimension), metadata(dat, 1))
    @test isnothing(metadata(dat, 1)[:window])
    @test label(dims(dat, 1)) == ""
    @test isnothing(units(dims(dat, 1)))
    @test metadata(dat, 1)[:test] == 123
    @test metadata(dims(dat, 1), :test) == 123

    @test metadatahelp(:window) == "Window function"
end

@testset "NMRBase: annotations" begin
    # Test with 19F CEST data
    cest_data = loadnmr(joinpath("test-data", "19f-cest-ts3"))

    # Test basic annotation access - returns entire annotations dictionary
    all_annot = annotations(cest_data)
    @test all_annot isa Dict{String,Any}
    @test haskey(all_annot, "title")

    # Test single key access with string
    @test annotations(cest_data, "title") == "19F CEST"

    # Test single key access with symbol
    @test annotations(cest_data, :title) == "19F CEST"

    # Test array access
    exp_type = annotations(cest_data, "experiment_type")
    @test exp_type isa Vector
    @test "cest" in exp_type

    # Test dimensions array
    dims_array = annotations(cest_data, "dimensions")
    @test dims_array[2] == "cest.offset"
    @test annotations(cest_data, "dimensions", 2) == "cest.offset"

    # Test nested dictionary access
    spinlock = annotations(cest_data, "cest")
    @test spinlock isa Dict
    @test annotations(cest_data, "cest", "channel") == "19F"
    @test annotations(cest_data, "cest", "duration") == 1
    @test annotations(cest_data, "cest.offset") isa FQList

    # Test FQList conversion functions
    fqlist = annotations(cest_data, "cest.offset")
    dim = dims(cest_data, F1Dim)

    # Test hz() for FQList (Hz, relative) - should return values as-is since they're already Hz relative
    hz_values = hz(fqlist, dim)
    @test hz_values[1] ≈ -2500.0  # First value
    @test hz_values[51] ≈ 0.0     # Middle value (on-resonance)
    @test hz_values[end] ≈ 2500.0 # Last value

    # Test ppm() for FQList - should convert to absolute ppm
    ppm_values = ppm(fqlist, dim)
    bf = dim[:bf]
    offsetppm = dim[:offsetppm]
    # Expected: offsetppm + Hz/bf*1e6
    @test ppm_values[1] ≈ offsetppm + (-2500.0 / bf * 1e6)
    @test ppm_values[51] ≈ offsetppm  # On-resonance point
    @test ppm_values[end] ≈ offsetppm + (2500.0 / bf * 1e6)

    # Test array of dictionaries (reference_pulse)
    refs = annotations(cest_data, "reference_pulse")
    @test refs isa Vector
    @test length(refs) == 1
    @test refs[1] isa Dict
    @test refs[1]["duration"] ≈ 13.29e-6
    @test annotations(cest_data, :reference_pulse, 1)["duration"] ≈ 13.29e-6
    p, pl = referencepulse(cest_data, "19F")
    @test p ≈ 13.29e-6
    @test watts(pl) ≈ 8.0

    # Test automatic dimension annotation - CEST data should get OffsetDim
    @test dims(cest_data, 2) isa OffsetDim
    @test cest_data[2, :units] == "ppm"
    @test cest_data[2, :nucleus] == "19F"
    # Verify offset values are in ppm (converted from FQList Hz)
    offset_axis = collect(dims(cest_data, 2))
    @test length(offset_axis) == 101  # 101 offset points
    # Values should be in ppm, not Hz (original was -2500 to +2500 Hz)
    @test all(abs.(offset_axis) .< 100)  # ppm values are much smaller than Hz

    # Test R1rho data - should get SpinlockDim and TrelaxDim
    r1rho_data = loadnmr(joinpath("test-data", "19f-r1rho-onres-ts4"))

    # Check dimensions array in annotations (should be reversed to match Julia order)
    dims_array = annotations(r1rho_data, "dimensions")
    @test dims_array[1] == "19F"
    @test dims_array[2] == "r1rho.power"
    @test dims_array[3] == "r1rho.duration"

    # Dimension 2 should be SpinlockDim with Hz values
    @test dims(r1rho_data, 2) isa SpinlockDim
    @test r1rho_data[2, :units] == "Hz"
    power_axis = collect(dims(r1rho_data, 2))
    @test length(power_axis) == 18  # 18 power levels
    @test all(power_axis .> 0)  # All positive Hz values
    @test power_axis[1] < power_axis[end]  # Increasing field strength

    # Dimension 3 should be TrelaxDim with seconds values
    @test dims(r1rho_data, 3) isa TrelaxDim
    @test r1rho_data[3, :units] == "s"
    @test r1rho_data[3, :delay_type] == :r1rho
    duration_axis = collect(dims(r1rho_data, 3))
    @test length(duration_axis) == 11  # 11 duration points
    @test duration_axis[1] ≈ 1e-4 atol = 1e-5  # First duration ~100 µs
    @test duration_axis[end] ≈ 0.08 atol = 0.01  # Last duration ~80 ms
end

@testset "NMRBase: Programmatic list annotations" begin
    # Create a mock NMRData object for testing
    using NMRTools.NMRIO: resolve_programmatic_lists!

    # Test 1: Linear spacing with start/step
    spec1 = exampledata("pseudo2D_XSTE")  # 2k x 10
    annotations1 = Dict{String,Any}("dimensions" => ["calibration.duration", "f1"],
                                    "calibration" =>
                                        Dict{String,Any}("duration" => Dict{String,
                                                                            Any}("type" => "linear",
                                                                                        "start" => 0.001,
                                                                                        "step" => 0.002)))
    resolve_programmatic_lists!(annotations1, spec1)

    # Should create a vector of length 206 (first dimension)
    @test annotations1["calibration"]["duration"] isa Vector
    @test length(annotations1["calibration"]["duration"]) == 10
    @test annotations1["calibration"]["duration"][1] ≈ 0.001
    @test annotations1["calibration"]["duration"][2] ≈ 0.003
    @test annotations1["calibration"]["duration"][end] ≈ 0.019

    # Test 2: Linear spacing with start/end
    spec2 = exampledata("pseudo2D_XSTE")  # 2k x 10
    annotations2 = Dict{String,Any}("dimensions" => ["calibration.duration", "f1"],
                                    "calibration" =>
                                        Dict{String,Any}("duration" => Dict{String,
                                                                            Any}("type" => "linear",
                                                                                        "start" => 0.001,
                                                                                        "end" => 0.1)))
    resolve_programmatic_lists!(annotations2, spec2)

    @test annotations2["calibration"]["duration"] isa Vector
    @test length(annotations2["calibration"]["duration"]) == 10
    @test annotations2["calibration"]["duration"][1] ≈ 0.001
    @test annotations2["calibration"]["duration"][end] ≈ 0.1
    # Check it's evenly spaced
    spacing = (0.1 - 0.001) / (10 - 1)
    @test annotations2["calibration"]["duration"][2] ≈ 0.001 + spacing

    # Test 3: Logarithmic spacing with start/end
    spec3 = exampledata("pseudo2D_XSTE")  # 2k x 10
    annotations3 = Dict{String,Any}("dimensions" => ["calibration.duration", "f1"],
                                    "calibration" =>
                                        Dict{String,Any}("duration" => Dict{String,
                                                                            Any}("type" => "log",
                                                                                        "start" => 0.001,
                                                                                        "end" => 0.1)))
    resolve_programmatic_lists!(annotations3, spec3)

    @test annotations3["calibration"]["duration"] isa Vector
    @test length(annotations3["calibration"]["duration"]) == 10
    @test annotations3["calibration"]["duration"][1] ≈ 0.001
    @test annotations3["calibration"]["duration"][end] ≈ 0.1
    # Check it's logarithmically spaced
    @test log(annotations3["calibration"]["duration"][2]) ≈
          log(0.001) + (log(0.1) - log(0.001)) / (10 - 1)

    # Test 4: Second dimension mapping
    spec4 = exampledata("pseudo2D_XSTE")  # 2k x 10
    annotations4 = Dict{String,Any}("dimensions" => ["gradient.strength", "f1"],
                                    "gradient" =>
                                        Dict{String,Any}("strength" => Dict{String,
                                                                            Any}("type" => "linear",
                                                                                        "start" => 0.05,
                                                                                        "step" => 0.1)))
    resolve_programmatic_lists!(annotations4, spec4)

    # Should create a vector of length 10 (second dimension)
    @test annotations4["gradient"]["strength"] isa Vector
    @test length(annotations4["gradient"]["strength"]) == 10
    @test annotations4["gradient"]["strength"][1] ≈ 0.05
    @test annotations4["gradient"]["strength"][2] ≈ 0.15
    @test annotations4["gradient"]["strength"][end] ≈ 0.95

    # Test 5: Counter-based pattern with integer counter
    spec5 = exampledata("pseudo2D_XSTE")
    annotations5 = Dict{String,Any}("dimensions" => ["calibration.duration", "f1"],
                                    "calibration" =>
                                        Dict{String,Any}("duration" => Dict{String,
                                                                            Any}("counter" => [0,
                                                                                               1,
                                                                                               2,
                                                                                               3,
                                                                                               4,
                                                                                               5,
                                                                                               6,
                                                                                               7,
                                                                                               8,
                                                                                               9],
                                                                                        "scale" => 0.01)))
    resolve_programmatic_lists!(annotations5, spec5)

    @test annotations5["calibration"]["duration"] isa Vector
    @test length(annotations5["calibration"]["duration"]) == 10
    @test annotations5["calibration"]["duration"][1] ≈ 0.0
    @test annotations5["calibration"]["duration"][2] ≈ 0.01
    @test annotations5["calibration"]["duration"][5] ≈ 0.04
    @test annotations5["calibration"]["duration"][end] ≈ 0.09

    # Test 6: Counter-based pattern with non-uniform counter
    spec6 = exampledata("pseudo2D_XSTE")
    annotations6 = Dict{String,Any}("dimensions" => ["gradient.strength", "f1"],
                                    "gradient" =>
                                        Dict{String,Any}("strength" => Dict{String,
                                                                            Any}("counter" => [1,
                                                                                               2,
                                                                                               4,
                                                                                               8,
                                                                                               16],
                                                                                        "scale" => 0.05)))
    resolve_programmatic_lists!(annotations6, spec6)

    @test annotations6["gradient"]["strength"] isa Vector
    @test length(annotations6["gradient"]["strength"]) == 5
    @test annotations6["gradient"]["strength"][1] ≈ 0.05
    @test annotations6["gradient"]["strength"][2] ≈ 0.10
    @test annotations6["gradient"]["strength"][3] ≈ 0.20
    @test annotations6["gradient"]["strength"][4] ≈ 0.40
    @test annotations6["gradient"]["strength"][5] ≈ 0.80
end

@testset "NMRBase: FQList conversions" begin
    # Create a test axis with known parameters
    spec = exampledata("1D_1H")
    ax = dims(spec, 1)
    bf = ax[:bf]
    offsetppm = ax[:offsetppm]

    # Test 1: Hz, relative (most common case)
    fq_hz_rel = FQList([-100.0, 0.0, 100.0], :Hz, true)
    @test hz(fq_hz_rel, ax) ≈ [-100.0, 0.0, 100.0]
    @test ppm(fq_hz_rel, ax) ≈ [offsetppm - 100.0 / bf * 1e6,
                                offsetppm,
                                offsetppm + 100.0 / bf * 1e6]

    # Test 2: Hz, absolute
    test_hz_abs = [bf * 7.0 * 1e-6,
                   bf * 7.5 * 1e-6,
                   bf * 8.0 * 1e-6]
    fq_hz_abs = FQList(test_hz_abs, :Hz, false)
    @test ppm(fq_hz_abs, ax) ≈ [7.0, 7.5, 8.0]
    @test hz(fq_hz_abs, ax) ≈ test_hz_abs .- ax[:offsethz]

    # Test 3: ppm, relative
    fq_ppm_rel = FQList([0.1, 0.0, -0.1], :ppm, true)
    @test ppm(fq_ppm_rel, ax) ≈ [offsetppm + 0.1, offsetppm, offsetppm - 0.1]
    @test hz(fq_ppm_rel, ax) ≈ [0.1 * bf * 1e-6, 0.0, -0.1 * bf * 1e-6]

    # Test 4: ppm, absolute (chemical shifts)
    fq_ppm_abs = FQList([7.0, 7.5, 8.0], :ppm, false)
    @test ppm(fq_ppm_abs, ax) ≈ [7.0, 7.5, 8.0]
    @test hz(fq_ppm_abs, ax) ≈ [(7.0 - offsetppm) * bf * 1e-6,
                                (7.5 - offsetppm) * bf * 1e-6,
                                (8.0 - offsetppm) * bf * 1e-6]

    # Test length
    @test length(fq_hz_rel) == 3
    @test length(fq_ppm_abs) == 3

    # Axis-free conversions using a reference baked into the list
    bf_ref = 6.0e8
    sfo_ref = bf_ref * (1 + 4.7e-6)   # carrier 4.7 ppm above bf
    offppm_ref = (sfo_ref - bf_ref) * 1e6 / bf_ref

    fq_free = FQList([-100.0, 0.0, 100.0], :Hz, true, bf_ref, sfo_ref)
    @test hz(fq_free) ≈ [-100.0, 0.0, 100.0]
    @test ppm(fq_free) ≈ [offppm_ref - 100.0 / bf_ref * 1e6,
                          offppm_ref,
                          offppm_ref + 100.0 / bf_ref * 1e6]

    # withreference attaches a reference from a channel dictionary
    ch = Dict{Symbol,Any}(:bf => bf_ref, :sfo => sfo_ref)
    fq_ref = withreference(FQList([-100.0, 0.0, 100.0], :Hz, true), ch)
    @test ppm(fq_ref) ≈ ppm(fq_free)
    @test hz(fq_ref) ≈ hz(fq_free)

    # A list with no reference cannot self-convert
    fq_noref = FQList([1.0, 2.0], :Hz, false)
    @test_throws NMRToolsError ppm(fq_noref)
    @test_throws NMRToolsError hz(fq_noref)
end

@testset "NMRBase: spectrometer channels" begin
    spec = exampledata("2D_HN")

    ch = channels(spec)
    @test ch isa AbstractDict
    # HN experiment routes 1H/13C/15N on f1/f2/f3 even though only H and N are detected
    @test Set(keys(ch)) == Set([:f1, :f2, :f3])

    # lookup by label
    @test channels(spec, :f1)[:nucleus] == H1
    @test channels(spec, :f3)[:nucleus] == N15
    @test channels(spec, :f9) === nothing

    # lookup by nucleus
    @test channels(spec, H1)[:label] == :f1
    @test channels(spec, N15)[:label] == :f3

    # carbon has no detected axis but is reachable through the channel model
    cch = channels(spec, C13)
    @test cch[:label] == :f2
    @test cch[:bf] == acqus(spec, :bf2)
    @test cch[:sfo] == acqus(spec, :sfo2)
    @test cch[:offsethz] ≈ cch[:sfo] - cch[:bf]
    @test cch[:offsetppm] ≈ (cch[:sfo] - cch[:bf]) * 1e6 / cch[:bf]

    # lookup by string (label or nucleus notation)
    @test channels(spec, "f2") === cch
    @test channels(spec, "13C") === cch
    @test channels(spec, "not-a-nucleus") === nothing
end

@testset "NMRBase: nuclei and coherences" begin
    @testset "parse nucleus tests" begin
        @testset "Mass number followed by element symbol" begin
            @test nucleus("1H") == H1
            @test nucleus("2H") == H2
            @test nucleus("12C") == C12
            @test nucleus("13C") == C13
            @test nucleus("14N") == N14
            @test nucleus("15N") == N15
            @test nucleus("19F") == F19
            @test nucleus("31P") == P31
        end

        @testset "Element symbol followed by mass number" begin
            @test nucleus("H1") == H1
            @test nucleus("H2") == H2
            @test nucleus("C12") == C12
            @test nucleus("C13") == C13
            @test nucleus("N14") == N14
            @test nucleus("N15") == N15
            @test nucleus("F19") == F19
            @test nucleus("P31") == P31
        end

        @testset "Case insensitive parsing" begin
            @test nucleus("19f") == F19
            @test nucleus("1h") == H1
            @test nucleus("13c") == C13
            @test nucleus("h1") == H1
            @test nucleus("c13") == C13
            @test nucleus("f19") == F19
        end

        @testset "Whitespace handling" begin
            @test nucleus(" 19F ") == F19
            @test nucleus("  1H  ") == H1
            @test nucleus("\t13C\n") == C13
        end

        @testset "Invalid nucleus strings" begin
            @test_throws ArgumentError nucleus("16O")  # Not defined in enum
            @test_throws ArgumentError nucleus("23Na") # Not defined in enum
            @test_throws ArgumentError nucleus("O16")  # Not defined in enum
            @test_throws ArgumentError nucleus("Na23") # Not defined in enum
        end

        @testset "Invalid format strings" begin
            @test_throws ArgumentError nucleus("H")     # Missing mass number
            @test_throws ArgumentError nucleus("13")    # Missing element
            @test_throws ArgumentError nucleus("ABC13") # Invalid element
            @test_throws ArgumentError nucleus("13XYZ") # Invalid element
            @test_throws ArgumentError nucleus("")      # Empty string
            @test_throws ArgumentError nucleus("H1C")   # Invalid format
            @test_throws ArgumentError nucleus("1H2")   # Invalid format
        end

        @testset "Edge cases" begin
            @test_throws ArgumentError nucleus("0H")    # Zero mass number
            @test_throws ArgumentError nucleus("H0")    # Zero mass number
            @test_throws ArgumentError nucleus("1000H") # Very large mass number
        end
    end # nucleus tests

    @test spin(H1) == 1 // 2
    @test spin(H2) == 1
    @test spin(C12) == 0
    @test gyromagneticratio(H1) == 2.6752218744e8

    @test coherenceorder(SQ(H1)) == 1
    @test coherenceorder(MQ(((H1, 1), (C13, 1)))) == 2
    @test coherenceorder(MQ(((H1, 1), (C13, -1)))) == 0
    @test coherenceorder(MQ(((H1, 3), (C13, 1)))) == 4
    @test coherenceorder(MQ(((H1, 0),))) == 0

    @test gyromagneticratio(SQ(H1)) == 2.6752218744e8
    @test gyromagneticratio(MQ(((H1, 1), (C13, 1)))) == 3.3480498744e8
    @test gyromagneticratio(MQ(((H1, 1), (C13, -1)))) == 2.0023938744e8
    @test gyromagneticratio(MQ(((H1, 3), (C13, 1)))) == 8.698493623199999e8
    @test gyromagneticratio(MQ(((H1, 0),))) == 0
end

@testset "NMRBase: window functions" begin
    dat = exampledata("pseudo2D_XSTE")
    @test dat[1, :window] isa Cos2Window
    @test dat[1, :window].tmax == 0.0512
    @test length(apod(dat, 1)) == 2048
    @test apod(dat, 1)[10] == 0.9992377902866475

    @test lineshape(dims(dat, 1), 8.0, 5)[1] ≈ 1.1210822807171034e-8
    @test lineshape(dims(dat, 1), 8.0, 5)[682] ≈ -0.00011946940995147479
    @test lineshape(dims(dat, 1), 8.0, 5)[686] ≈ 0.02244850216794408
    @test lineshape(dims(dat, 1), 8.0, 5, ComplexLineshape())[1] ≈
          1.1210822807171034e-8 + 4.754164488338695e-5im
end

@testset "NMRBase: decimate" begin
    # generic function (not NMR data)
    @test decimate(1:10, 2) == [1.5, 3.5, 5.5, 7.5, 9.5]

    # NMR data
    dat = exampledata("pseudo2D_XSTE")
    @test size(decimate(dat, 10, 1)) == (204, 10)
    @test decimate(dat, 10, 1)[2, 2] == 43.59735595703125
    @test size(decimate(dat, 3, 2)) == (2048, 3)
    @test decimate(dat, 3, 2)[2, 2] == -31.345370822482636
end

@testset "NMRBase: stack" begin
    specs = [loadnmr(joinpath(artifact"1D_19F_titration", "1")),
             loadnmr(joinpath(artifact"1D_19F_titration", "2")),
             loadnmr(joinpath(artifact"1D_19F_titration", "3"))]
    stacked = stack(specs)
    @test size(stacked) == (32768, 3)

    specs2 = [loadnmr(joinpath(artifact"1D_19F_titration", "1")),
              exampledata("2D_HN")]
    @test_throws DimensionMismatch stack(specs2)

    specs3 = [loadnmr(joinpath(artifact"2D_HN_titration", "1/test.ft2")),
              loadnmr(joinpath(artifact"2D_HN_titration", "2/test.ft2")),
              loadnmr(joinpath(artifact"2D_HN_titration", "3/test.ft2")),
              loadnmr(joinpath(artifact"2D_HN_titration", "4/test.ft2"))]
    # expect warning that experiments do not have same rg
    @test (@test_logs (:warn,) size(stack(specs3))) == (768, 512, 4)
end

@testset "NMRBase: Power" begin
    # Test construction from dB
    p1 = Power(30.0, :dB)
    @test db(p1) == 30.0
    @test watts(p1) ≈ 0.001

    # Test construction from Watts
    p2 = Power(1.0, :W)
    @test watts(p2) == 1.0
    @test db(p2) == 0.0

    # Test with integer input
    p3 = Power(10, :W)
    @test watts(p3) == 10.0
    @test db(p3) ≈ -10.0

    # Test zero Watts handling
    p4 = Power(0.0, :W)
    @test db(p4) == 120.0
    @test watts(p4) ≈ 1e-12

    # Test error for missing unit
    @test_throws ArgumentError Power(10)

    # Test error for invalid unit
    @test_throws ArgumentError Power(10, :V)

    # Test roundtrip conversions
    original_watts = 0.5
    p5 = Power(original_watts, :W)
    @test watts(p5) ≈ original_watts

    original_db = -6.0
    p6 = Power(original_db, :dB)
    @test db(p6) == original_db

    # Test pretty printing contains both units
    p7 = Power(20.0, :dB)
    output = string(p7)
    @test occursin("dB", output)
    @test occursin("W", output)

    # Test hz(::Power, ref_p, ref_Hz) - same power returns same Hz
    ref = Power(10.0, :dB)
    @test hz(ref, ref, 500.0) ≈ 500.0

    # 6 dB less attenuation → 2x the field strength
    p8 = Power(4.0, :dB)
    @test hz(p8, ref, 500.0) ≈ 500.0 * 10^(6.0 / 20) atol = 0.01

    # Test hz(::Power, ref_p, ref_pulselength, ref_angle) variant
    # 10 μs 90° pulse → 1/(4*10e-6) = 25000 Hz
    ref_pl = Power(10.0, :dB)
    @test hz(ref_pl, ref_pl, 10e-6, 90) ≈ 25000.0

    # Consistency: both hz() methods should agree
    ref_hz_direct = hz(p8, ref, 25000.0)
    ref_hz_pulse = hz(p8, ref, 10e-6, 90)
    @test ref_hz_direct ≈ ref_hz_pulse
end

@testset "NMRBase: Peak type" begin
    p = Peak(8.2, 120.5; label="G10")
    @test p isa Peak{2}
    @test ppm(p) == (8.2, 120.5)
    @test ppm(p, 2) == 120.5
    @test label(p) == "G10"
    @test isnan(intensity(p))
    @test all(isnan, linewidth(p))
    @test ndims(p) == 2

    q = Peak(p; label="G10N", intensity=2.0)
    @test label(q) == "G10N" && intensity(q) == 2.0 && ppm(q) == ppm(p)
    @test Peak(p) == p
    @test Peak((4.7,); linewidth=(1.5,)) == Peak(4.7; linewidth=(1.5,))
    @test_throws DimensionMismatch Peak(8.2, 120.5; linewidth=(1.0,))

    pks = [Peak(8.2, 120.5; label="A"), Peak(7.1, 115.0; intensity=3.0)]
    cols = NMRTools.NMRBase.Tables.columns(pks)
    @test cols.label == ["A", ""]
    @test cols.x == [8.2, 7.1]
    @test cols.y == [120.5, 115.0]
    @test isequal(cols.intensity, [NaN, 3.0])
    @test haskey(cols, :lwx) && haskey(cols, :lwy)
    @test occursin("2-element Vector{Peak{2}}", repr("text/plain", pks))
end

@testset "NMRBase: findnearest" begin
    pks = [Peak(8.0, 120.0), Peak(8.5, 110.0), Peak(7.0, 121.0)]
    @test findnearest(pks, (8.4, 119.0)) == 1
    # with ¹⁵N down-weighted, the ¹H shift dominates
    @test findnearest(pks, (8.4, 112.0); scale=(1, 0.01)) == 2
    @test findnearest(pks, (7.2, 116.0); scale=(1, 0.01)) == 3
    @test isnothing(findnearest(Peak{2}[], (8.0, 120.0)))
    @test findnearest([Peak(1.0), Peak(2.0)], 1.8) == 2
    @test_throws DimensionMismatch findnearest(pks, 8.0)
end

@testset "NMRBase: Peak detection (1D)" begin
    spec = exampledata("1D_19F")
    pks = detectpeaks(spec)
    @test pks isa Vector{Peak{1}}
    @test length(pks) == 1
    @test label(pks[1]) == "X1"
    @test ppm(pks[1], 1) ≈ -123.698 atol = 0.005
    @test intensity(pks[1]) > 0
    @test 1 < linewidth(pks[1], 1) < 100   # Hz
    @test length(detectpeaks(spec; snr=0.1)) > length(pks)
    @test isempty(detectpeaks(spec; snr=1e6))
    # deprecated keyword from v0.2
    @test isempty(@test_deprecated detectpeaks(spec; snr_threshold=1e6))
    @test isempty(detectpeaks(spec; sign=:negative))
    @test_throws ArgumentError detectpeaks(spec; sign=:up)

    # an absolute threshold is the same as snr × noise
    @test detectpeaks(spec; threshold=10 * spec[:noise]) == pks

    # findpeak climbs to the same maximum from nearby
    p = findpeak(spec, ppm(pks[1], 1) + 0.02; label="F")
    @test ppm(p, 1) ≈ ppm(pks[1], 1)
    @test label(p) == "F"

    # pseudo-2D data must be sliced first
    @test_throws ArgumentError detectpeaks(exampledata("pseudo2D_XSTE"))
end

@testset "NMRBase: Peak detection (synthetic Lorentzian)" begin
    spec = exampledata("1D_19F")
    ax = dims(spec, F1Dim)
    ν = hz(ax)
    for (δ0, fwhm) in ((-123.7, 20.0), (-121.0, 50.0))
        ν0 = hz(δ0, ax)
        y = @. 1e6 / (1 + ((ν - ν0) / (fwhm / 2))^2)
        synth = NMRData(y, dims(spec); metadata=deepcopy(metadata(spec)))
        synth[:noise] = 1.0
        pks = detectpeaks(synth)
        @test length(pks) == 1
        @test ppm(pks[1], 1) ≈ δ0 atol = 1e-4
        @test intensity(pks[1]) ≈ 1e6 rtol = 0.02
        @test linewidth(pks[1], 1) ≈ fwhm rtol = 0.05
        @test findpeak(synth, δ0 + 0.1) == Peak(pks[1]; label="")

        neg = NMRData(-y, dims(spec); metadata=deepcopy(metadata(spec)))
        neg[:noise] = 1.0
        @test intensity(only(detectpeaks(neg))) ≈ -1e6 rtol = 0.02
        @test isempty(detectpeaks(neg; sign=:positive))
    end
end

@testset "NMRBase: Peak detection (2D)" begin
    spec = exampledata("2D_HN")
    pks = detectpeaks(spec; sign=:positive)
    @test pks isa Vector{Peak{2}}
    @test length(pks) > 10
    @test issorted(intensity.(pks); rev=true)
    @test all(p -> intensity(p) > 10 * spec[:noise], pks)
    xr, yr = extrema(data(dims(spec, 1))), extrema(data(dims(spec, 2)))
    @test all(p -> xr[1] ≤ ppm(p, 1) ≤ xr[2] && yr[1] ≤ ppm(p, 2) ≤ yr[2], pks)
    @test all(p -> all(lw -> isnan(lw) || 0 < lw < 500, linewidth(p)), pks)
    @test length(unique(label.(pks))) == length(pks)

    # findpeak from slightly off the strongest peak returns it
    p1 = pks[1]
    p = findpeak(spec, ppm(p1) .+ (0.01, 0.1))
    @test ppm(p) == ppm(p1)
    @test intensity(p) == intensity(p1)
    @test_throws DimensionMismatch findpeak(spec, (8.0,))
end

@testset "NMRBase: register!" begin
    spec = exampledata("2D_HN")
    step = abs.((data(dims(spec, 1))[2] - data(dims(spec, 1))[1],
                 data(dims(spec, 2))[2] - data(dims(spec, 2))[1]))
    pks = detectpeaks(spec; sign=:positive)[1:40]
    offset = (0.12, -1.3)
    moved = [Peak(p; position=ppm(p) .+ offset) for p in pks]

    # the non-mutating form leaves the list alone
    back = register(moved, spec)
    @test ppm(moved[1]) == ppm(pks[1]) .+ offset
    @test all(i -> all(abs.(ppm(back[i]) .- ppm(pks[i])) .≤ step), eachindex(pks))

    # register! moves the list in place and returns it
    @test register!(moved, spec) === moved
    @test moved == back

    # a dimension can be held fixed
    held = register([Peak(p; position=ppm(p) .+ (0.12, 0.0)) for p in pks], spec;
                    maxshift=(0.3, 0))
    @test all(i -> ppm(held[i], 2) == ppm(pks[i], 2), eachindex(pks))
    @test all(i -> abs(ppm(held[i], 1) - ppm(pks[i], 1)) ≤ step[1], eachindex(pks))

    # an offset beyond maxshift cannot be found
    far = register([Peak(p; position=ppm(p) .+ (0.12, 0.0)) for p in pks], spec;
                   maxshift=0.05)
    @test abs(ppm(far[1], 1) - ppm(pks[1], 1)) > 0.05

    # negative peaks count as much as positive ones, and one huge peak elsewhere
    # does not outweigh several weaker matches
    spec1 = exampledata("1D_19F")
    ax = dims(spec1, F1Dim)
    ν = hz(ax)
    function lorentz(δ, h)
        ν0 = hz(δ, ax)
        return @. h / (1 + ((ν - ν0) / 10)^2)
    end
    # a peak 100 times stronger, which an offset of +0.6 ppm would land the first peak on
    y = lorentz(-120.5, -1e4) .+ lorentz(-121.5, -1e4) .+ lorentz(-122.5, -1e4) .+
        lorentz(-119.7, 1e6)
    synth = NMRData(y, dims(spec1); metadata=deepcopy(metadata(spec1)))
    synth[:noise] = 10.0
    predicted = [Peak(δ + 0.2) for δ in (-120.5, -121.5, -122.5)]
    register!(predicted, synth; maxshift=1.0)
    @test all(isapprox.(ppm.(predicted, 1), [-120.5, -121.5, -122.5]; atol=0.005))

    @test isempty(register!(Peak{2}[], spec))
    @test_throws DimensionMismatch register(pks, spec; maxshift=(0.1, 0.1, 0.1))
end

@testset "NMRBase: Metadata independence after slicing (issue #25)" begin
    # Test that sliced NMRData has independent metadata from the original
    # See GitHub issue #25: https://github.com/waudbylab/NMRTools.jl/issues/25

    # 2D pseudo-data test
    dat1a = exampledata("pseudo2D_XSTE")
    dat3 = dat1a[:, 3]  # slice creates new NMRData

    # Set label on the slice
    label!(dat3, "sliced data")
    # Original should NOT be affected
    @test label(dat1a) != "sliced data"

    # Set label on the original
    label!(dat1a, "original data")
    # Slice should NOT be affected
    @test label(dat3) == "sliced data"

    # Test metadata modification via setindex!
    dat3[:noise] = 123.0
    @test dat1a[:noise] != 123.0  # original should not change

    # Test that NMRData constructor also creates independent copy
    dat2 = exampledata("1D_19F")
    dat2b = NMRData(dat2)
    label!(dat2b, "copy label")
    @test label(dat2) != "copy label"  # original should not change

    # Test slicing via interval selector
    dat1c = dat1a[7.0 .. 9.0, :]
    label!(dat1c, "interval slice")
    @test label(dat1a) != "interval slice"

    # Test selectdim behavior (used by decimate)
    dat1d = selectdim(dat1a, 2, 1:5)
    label!(dat1d, "selectdim slice")
    @test label(dat1a) != "selectdim slice"

    # Test estimatenoise doesn't affect original after slicing
    dat_for_noise = exampledata("1D_1H")
    dat_slice = dat_for_noise[7.0 .. 9.0]
    estimatenoise(dat_slice)
    # The original should have different noise value (or nil if not yet estimated)
    @test dat_for_noise[:noise] != dat_slice[:noise]
end

@testset "NMRBase: Xi ratios" begin
    # Test Xi ratio retrieval for DSS (aqueous)
    @test xi_ratio(H1) == 1.0
    @test xi_ratio(H1; aqueous=true) == 1.0
    @test xi_ratio(C13) == 0.251449530
    @test xi_ratio(C13; aqueous=true) == 0.251449530

    # Test Xi ratio retrieval for TMS (organic)
    @test xi_ratio(H1; aqueous=false) == 1.0
    @test xi_ratio(C13; aqueous=false) == 0.25145020

    # Test undefined nucleus returns nothing
    @test xi_ratio(C12) === nothing
    @test xi_ratio(N14) === nothing
end

@testset "NMRBase: finddim and resolvedim" begin
    spec = exampledata("1D_1H")

    # Test finddim with Nucleus
    @test finddim(spec, H1) == 1
    @test finddim(spec, C13) === nothing  # No C13 dimension

    # Test resolvedim with different input types
    @test resolvedim(spec, 1) == 1
    @test resolvedim(spec, F1Dim) == 1
    @test resolvedim(spec, H1) == 1

    # Test resolvedim error for non-existent nucleus
    @test_throws ArgumentError resolvedim(spec, C13)
end

@testset "NMRBase: shiftdim with nucleus" begin
    spec = exampledata("1D_1H")
    original_ppm = data(spec, 1)[1]

    # Test shiftdim with nucleus
    spec2 = shiftdim(spec, H1, 0.5)
    @test data(spec2, 1)[1] - original_ppm ≈ 0.5

    # Verify referenceoffset metadata is updated
    @test metadata(spec2, 1, :referenceoffset) ≈ 0.5
end

@testset "NMRBase: watershift" begin
    # Test water shift calculation at various temperatures
    # Formula: δ(H2O) = 7.83 − T / 96.9

    # At 300 K (approx 27°C)
    @test watershift(300.0) ≈ 7.83 - 300.0 / 96.9

    # At 273.15 K (0°C)
    @test watershift(273.15) ≈ 7.83 - 273.15 / 96.9

    # At 298.15 K (25°C)
    @test watershift(298.15) ≈ 7.83 - 298.15 / 96.9

    # Verify the shift is around 4.7 ppm at physiological temperature
    @test watershift(310.0) ≈ 4.630958 atol = 0.01
end

@testset "NMRBase: isaqueous" begin
    spec = exampledata("1D_1H")
    @test isaqueous(spec) == true
end

@testset "NMRBase: reference function" begin
    spec = exampledata("1D_1H")
    original_ppm = data(spec, 1)[1]

    # Test basic referencing with explicit shift pair
    @test_logs (:info,) begin
        spec2 = reference(spec, 1, 4.7 => 4.8; indirect=false)
        @test data(spec2, 1)[1] - original_ppm ≈ 0.1
    end

    # Test referencing with F1Dim
    @test_logs (:info,) begin
        spec3 = reference(spec, F1Dim, 4.7 => 4.8; indirect=false)
        @test data(spec3, 1)[1] - original_ppm ≈ 0.1
    end

    # Test referencing with nucleus
    @test_logs (:info,) begin
        spec4 = reference(spec, H1, 4.7 => 4.8; indirect=false)
        @test data(spec4, 1)[1] - original_ppm ≈ 0.1
    end

    # Test referenceoffset is updated correctly
    @test_logs (:info,) begin
        spec6 = reference(spec, 1, 0.0 => 0.5; indirect=false)
        @test metadata(spec6, 1, :referenceoffset) ≈ 0.5
    end
end

@testset "NMRBase: reference with 2D data" begin
    spec = exampledata("2D_HN")

    # Get original ppm values
    original_f1 = data(spec, 1)[1]
    original_f2 = data(spec, 2)[1]

    # Test referencing single dimension without indirect
    @test_logs (:info,) begin
        spec2 = reference(spec, 1, 8.0 => 8.1; indirect=false)
        @test data(spec2, 1)[1] - original_f1 ≈ 0.1
        @test data(spec2, 2)[1] == original_f2  # F2 unchanged
    end

    # Test multi-dimension referencing
    @test_logs (:info,) (:info,) begin
        spec3 = reference(spec, [F1Dim, F2Dim], [8.0 => 8.1, 120.0 => 119.5];
                          indirect=false)
        @test data(spec3, 1)[1] - original_f1 ≈ 0.1
        @test data(spec3, 2)[1] - original_f2 ≈ -0.5
    end

    # Test multi-dimension with tuples
    @test_logs (:info,) (:info,) begin
        spec4 = reference(spec, (F1Dim, F2Dim), (8.0 => 8.1, 120.0 => 119.5);
                          indirect=false)
        @test data(spec4, 1)[1] - original_f1 ≈ 0.1
        @test data(spec4, 2)[1] - original_f2 ≈ -0.5
    end
end