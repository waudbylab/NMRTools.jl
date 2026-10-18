using NMRTools
using Artifacts
using LazyArtifacts
using Test
using Dates

@testset "NMRIO: 1D 19F (nmrPipe)" begin
    dat = loadnmr(joinpath(artifact"1D_19F", "test.ft1"))

    @test size(dat) == (8751,)
    @test length(dat) == 8751

    @test parent(dat) isa Vector{Float64}

    @test dat[1] == -3741.875
    @test dat[end] == 5973.625
    @test maximum(dat) == 1.274037375e6

    @test label(dat) == "Example 19F 1D spectrum"
    @test dat[:title] == """Example 19F 1D spectrum
    Sample details in second line
    and third line"""
    @test dat[:ns] == 64
    @test dat[:rg] == 76.02
    @test dat[:format] == :nmrpipe
    @test dat[:ndim] == 1

    @test scale(dat) == 4865.28

    @test acqus(dat, :pulprog) == "zg"
    @test acqus(dat, :p, 1) == 13.25e-6
    @test acqus(dat, "P", 2) == 22e-6
    @test acqus(dat, :d, 1) == 4
    @test acqus(dat)[:te] == 318.1509

    @test dims(dat, 1) isa FrequencyDimension
    @test length(dims(dat, 1)) == 8751
    @test dims(dat, 1)[1] == -119.99959383646768
    @test dims(dat, 1)[2] == -120.00050820180846
    @test dims(dat, 1)[end] == -128.00029056823928

    @test metadata(dat, 1, :sf) == 470.4608353134155e6
    @test dat[1, :bf] == 470.5220031738281e6
    @test label(dat, 1) == "19F"
    @test dat[1, :td] == 14097
    @test dat[1, :tdzf] == 32768
    @test dat[1, :region] == 5448:14198
    @test dat[1, :offsetppm] == -130
    @test dat[1, :window] == ExponentialWindow(5.0, 0.9999472073130639)
    @test dat[1, :swhz] == 14097.744257799108
    @test dat[1, :swppm] == 29.961923486479087
    @test dat[1, :pseudodim] == false
end

@testset "NMRIO: 2D HN (nmrPipe)" begin
    dat = loadnmr(joinpath(artifact"2D_HN", "test.ft2"))

    @test size(dat) == (568, 256)
    @test length(dat) == 145408

    @test parent(dat) isa Matrix{Float64}

    @test dat[1, 100] == 59228.95703125
    @test dat[20, 100] == 27169.9140625
    @test dat[1] == -48482.26171875
    @test dat[1, end] == 36429.125
    @test dat[end, end] == -23958.6171875
    @test dat[end] == -23958.6171875
    @test maximum(dat) == 7.984743e6

    @test dat[Near(8.3), Near(120)] == 81722.171875
    @test dat[Between(8.3, 8.35), Near(120)] == [86881.734375
                                                 115982.59375
                                                 115535.171875
                                                 84537.9296875
                                                 51392.3046875
                                                 81722.171875]
    @test dat[Near(8.4), Between(110, 111)] == [-22987.546875
                                                18746.140625
                                                35149.3359375
                                                21203.44140625
                                                -2393.46875
                                                -22528.533203125
                                                -27708.4296875
                                                -3501.19873046875
                                                39905.0
                                                63979.9765625
                                                51346.94140625
                                                33063.59375]
    @test dat[Between(8.6, 8.62), Between(108, 108.4)] ==
          [57461.02734375 25821.890625 -20170.22265625 -16362.2734375
           27847.72265625 9394.47265625 -34206.65625 -59513.51171875
           -34801.6171875 -85781.984375 -96177.78125 -83699.1171875]

    @test label(dat) == "13C,15N ubiquitin"
    @test dat[:title] == """13C,15N ubiquitin
    500 uM in 10% D2O, 20 mM Pi pH 6.5, 277 K, SOFAST-HMQC"""
    @test dat[:ns] == 2
    @test dat[:rg] == 203
    @test dat[:format] == :nmrpipe
    @test dat[:ndim] == 2

    @test scale(dat) == 406

    @test acqus(dat, :pulprog) == "sfhmqcf3gpph.cw"
    @test acqus(dat, :p, 1) ≈ 9.2e-6
    @test acqus(dat, "P", 2) ≈ 18.4e-6
    @test acqus(dat, :d, 1) == 0.1
    @test acqus(dat)[:te] == 276.9988

    @test dims(dat, 1) isa FrequencyDimension
    @test length(dims(dat, 1)) == 568
    @test dims(dat, 1)[1] == 10.996211600312328
    @test dims(dat, 1)[2] == 10.987407197534884
    @test dims(dat, 1)[end] == 6.004115225500486

    @test dat[1, :bf] == 600.2030029296875e6
    @test metadata(dat, 1, :sf) == 600.2059883394103e6
    @test label(dat, 1) == "HN"
    @test dat[1, :td] == 1024
    @test dat[1, :tdzf] == 2048
    @test dat[1, :region] == 341:908
    @test dat[1, :offsetppm] == 4.973999977111816
    @test dat[1, :window] == CosWindow(0.09461760226549189)
    @test dat[1, :swhz] == 10822.510563380281
    @test dat[1, :swppm] == 18.031416888209264
    @test dat[1, :pseudodim] == false

    @test dat[2, :bf] == 60.82500076293945e6
    @test metadata(dat, 2, :sf) == 60.832195448377445e6
    @test label(dat, 2) == "15N"
    @test dat[2, :td] == 128
    @test dat[2, :tdzf] == 256
    @test ismissing(dat[2, :region])
    @test dat[2, :offsetppm] == 118.28500366210938
    @test dat[2, :window] == CosWindow(0.09566717062018215)
    @test dat[2, :swhz] == 1337.9720458984375
    @test dat[2, :swppm] == 21.997074050406937
    @test dat[2, :pseudodim] == false

    @test dat[:noise] == 42605.50537796347
end

@testset "NMRIO: pseudo-2D XSTE (nmrPipe)" begin
    dat = loadnmr(joinpath(artifact"pseudo2D_XSTE", "test.ft1"))

    @test size(dat) == (206, 10)
    @test length(dat) == 2060

    @test parent(dat) isa Matrix{Float64}

    @test dat[100, 5] == 193505.375
    @test dat[1, end] == 413.376953125
    @test dat[1] == 8857.43359375
    @test dat[end] == -5065.54248046875
    @test maximum(dat) == 1.277118e6

    @test dat[Near(8.3), 1] == 635359.375
    @test dat[Between(8.3, 8.35), 9] == [241633.140625
                                         219759.296875]
    @test dat[Near(8.4), :] == [1.277118e6
                                1.247608625e6
                                1.1866e6
                                1.094307375e6
                                996278.0625
                                887989.875
                                775422.0
                                641461.4375
                                520499.28125
                                399827.28125]
    @test dat[Between(7.8, 7.86), 3:5] == [27620.0078125 14238.609375 15452.2421875
                                           15589.8359375 -10553.5 -1367.3984375
                                           14052.1484375 -4170.84375 -2667.4921875]

    @test label(dat) == "15N aSyn, 283 K - XSTE"
    @test dat[:title] == """15N aSyn, 283 K - XSTE
    DELTA = 100ms, delta = 4ms, G = 5-95%"""
    @test dat[:ns] == 16
    @test dat[:rg] == 128
    @test dat[:format] == :nmrpipe
    @test dat[:ndim] == 2

    @test scale(dat) == 2048

    @test acqus(dat, :pulprog) == "stebpgp1s19xn.jk"
    @test acqus(dat, :p, 1) == 9.25e-6
    @test acqus(dat, "P", 2) == 18.5e-6
    @test acqus(dat, :d, 1) == 1
    @test acqus(dat)[:te] == 283.0556

    @test dims(dat, 1) isa FrequencyDimension
    @test length(dims(dat, 1)) == 206
    @test dims(dat, 1)[1] == 9.995628344572216
    @test dims(dat, 1)[2] == 9.976091311209679
    @test dims(dat, 1)[end] == 5.990536505252104

    @test metadata(dat, 1, :sf) == 499.85444701842004e6
    @test dat[1, :bf] == 499.85198974609375e6
    @test label(dat, 1) == "1H"
    @test dat[1, :td] == 512
    @test dat[1, :tdzf] == 1024
    @test dat[1, :region] == 253:458
    @test dat[1, :offsetppm] == 4.915999889373779
    @test dat[1, :window] == CosWindow(0.0512)
    @test dat[1, :swhz] == 10000.0
    @test dat[1, :swppm] == 20.00592216323802
    @test dat[1, :pseudodim] == false

    @test metadata(dat, 2, :sf) === nothing
    @test dat[2, :bf] === nothing
    @test label(dat, 2) == "TAU"
    @test dat[2, :npoints] == 10
    @test dat[2, :td] === nothing
    @test dat[2, :tdzf] === nothing
    @test dat[2, :region] === nothing
    @test dat[2, :pseudodim] == true
end

@testset "NMRIO: 1D 19F (Bruker pdata)" begin
    dat = exampledata("1D_19F")

    @test size(dat) == (5000,)
    @test length(dat) == 5000

    @test parent(dat) isa Vector{Float64}

    @test dat[1] == 4671.5184326171875
    @test dat[end] == 5857.3370361328125
    @test maximum(dat) == 23737.674377441406

    @test label(dat) == "Example 19F 1D spectrum"
    @test dat[:title] == """Example 19F 1D spectrum
    Sample details in second line
    and third line"""
    @test dat[:ns] == 64
    @test dat[:rg] == 76.02
    @test dat[:format] == :pdata
    @test dat[:ndim] == 1

    @test scale(dat) == 4865.28

    @test acqus(dat, :pulprog) == "zg"
    @test acqus(dat, :p, 1) == 13.25e-6
    @test acqus(dat, "P", 2) == 22e-6
    @test acqus(dat, :d, 1) == 4
    @test acqus(dat)[:te] == 318.1509

    @test dims(dat, 1) isa FrequencyDimension
    @test length(dims(dat, 1)) == 5000
    @test dims(dat, 1)[1] == -119.5922
    @test dims(dat, 1)[2] == -119.59311424688988
    @test dims(dat, 1)[end] == -124.16252020246354

    @test metadata(dat, 1, :sf) == 470.5217922204685e6
    @test dat[1, :bf] == 470.582968e6
    @test label(dat, 1) == "19F"
    @test dat[1, :td] == 14097
    @test dat[1, :tdzf] == 32768
    @test dat[1, :region] == 5000:9999
    @test dat[1, :offsetppm] == -129.99998659428783
    @test dat[1, :window] == ExponentialWindow(10.0, 0.9999471999999971)
    @test dat[1, :swhz] == 14097.744360902296
    @test dat[1, :swppm] == 29.95804208728246
    @test dat[1, :pseudodim] == false

    # load complex spectrum (imaginary part conjugated relative to TopSpin's 1i)
    dat2 = loadnmr(artifact"1D_19F"; allcomponents=true)
    @test dat2[100] == Multicomplex(4920.725646972656, 17351.53515625)
    @test metadata(dat2, 1, :mcindex) == 1
    # loaded data can be phased like processed data
    @test realest(phase(phase(dat2, 1, 30.0), 1, -30.0)) ≈ realest(dat2)
end

@testset "NMRIO: 2D HN (Bruker pdata)" begin
    dat = exampledata("2D_HN")

    @test size(dat) == (1216, 512)
    @test length(dat) == 622592

    @test parent(dat) isa Matrix{Float64}

    @test dat[1, 100] == 587.3193359375
    @test dat[20, 100] == 277.4541015625
    @test dat[1] == 20.85546875
    @test dat[1, end] == 40.3671875
    @test dat[end, end] == -45.9990234375
    @test dat[end] == -45.9990234375
    @test maximum(dat) == 88924.693359375

    @test dat[Near(8.3), Near(120)] == 3525.2021484375
    @test dat[Between(8.3, 8.35), Near(120)] == [2630.798828125
                                                 2762.41015625
                                                 3282.8671875
                                                 4318.3232421875
                                                 5236.74609375
                                                 5653.0947265625
                                                 5525.8525390625
                                                 5008.009765625
                                                 4389.3798828125
                                                 3896.55859375
                                                 3525.2021484375]
    @test dat[Near(8.4), Between(110, 111)] == [922.9013671875
                                                -1304.8193359375
                                                -2347.005859375
                                                -2250.8642578125
                                                -1314.556640625
                                                20.0224609375
                                                1305.328125
                                                2209.80859375
                                                2585.904296875
                                                2469.98828125
                                                2025.5947265625
                                                1459.6162109375
                                                945.8330078125
                                                580.998046875
                                                381.263671875
                                                309.458984375
                                                313.517578125
                                                356.509765625
                                                427.232421875
                                                531.966796875
                                                676.92578125
                                                853.533203125
                                                1034.7265625
                                                1183.2421875]
    @test dat[Between(8.6, 8.62), Between(108, 108.4)] ==
          [1654.3662109375 1509.92578125 1408.123046875 1357.8310546875 1359.12890625 1402.8837890625 1471.8896484375 1543.5458984375 1593.7431640625
           1897.447265625 1729.6943359375 1606.208984375 1536.0478515625 1519.2578125 1546.73828125 1601.5107421875 1661.31640625 1702.26171875
           2452.357421875 2218.82421875 2009.0166015625 1839.666015625 1721.0146484375 1655.0849609375 1634.916015625 1645.1904296875 1664.50390625
           3220.08203125 2852.16015625 2480.65234375 2141.1640625 1864.3203125 1670.0283203125 1563.08984375 1531.59765625 1549.009765625
           4082.63671875 3541.2294921875 2965.853515625 2417.021484375 1950.8349609375 1607.763671875 1403.9716796875 1327.916015625 1343.6220703125]

    @test label(dat) == "13C,15N ubiquitin"
    @test dat[:title] == """13C,15N ubiquitin
    500 uM in 10% D2O, 20 mM Pi pH 6.5, 277 K, SOFAST-HMQC"""
    @test dat[:ns] == 2
    @test dat[:rg] == 203
    @test dat[:format] == :pdata
    @test dat[:ndim] == 2

    @test scale(dat) == 406

    @test acqus(dat, :pulprog) == "sfhmqcf3gpph.cw"
    @test acqus(dat, :p, 1) ≈ 9.2e-6
    @test acqus(dat, "P", 2) ≈ 18.4e-6
    @test acqus(dat, :d, 1) == 0.1
    @test acqus(dat)[:te] == 276.9988

    @test dims(dat, 1) isa FrequencyDimension
    @test length(dims(dat, 1)) == 1216
    @test dims(dat, 1)[1] == 11.07119
    @test dims(dat, 1)[2] == 11.066787776480712
    @test dims(dat, 1)[end] == 5.722488424064256

    @test metadata(dat, 1, :sf) == 600.2028190015606e6
    @test dat[1, :bf] == 600.2e6
    @test label(dat, 1) == "1H"
    @test dat[1, :td] == 1024
    @test dat[1, :tdzf] == 4096
    @test dat[1, :region] == 600:1815
    @test dat[1, :offsetppm] ≈ 4.696770344069996
    @test dat[1, :window] == Cos2Window(0.09461759999999972)
    @test dat[1, :swhz] == 10822.510822510854
    @test dat[1, :swppm] ≈ 18.03150753500642
    @test dat[1, :pseudodim] == false

    @test metadata(dat, 2, :sf) == 60.82491449022095e6
    @test dat[2, :bf] == 60.817738e6
    @test label(dat, 2) == "15N"
    @test dat[2, :td] == 64
    @test dat[2, :tdzf] == 512
    @test ismissing(dat[2, :region])
    @test dat[2, :offsetppm] == 117.99995292412437
    @test dat[2, :window] == CosWindow(0.04783359999999993)
    @test dat[2, :swhz] == 1337.97163500134
    @test dat[2, :swppm] == 21.999694151751253
    @test dat[2, :pseudodim] == false

    # field strength metadata (in Tesla)
    @test dat[:field] isa Float64
    @test 14.0 < dat[:field] < 14.2  # ~600 MHz spectrometer ≈ 14.1 T

    # load multicomplex spectrum (bicomplex for 2D)
    datMC = loadnmr(artifact"2D_HN"; allcomponents=true)
    @test size(datMC) == (1216, 512)
    @test parent(datMC) isa Matrix{<:Multicomplex}
    @test order(datMC[1, 1]) == 2  # bicomplex
    @test realest(datMC[1, 100]) == dat[1, 100]  # real-real component matches real-only import
    @test metadata(datMC, 1, :mcindex) == 1 && metadata(datMC, 2, :mcindex) == 2
    @test metadata(dat, 2, :mcindex) == 2         # set for real-only data too
end

@testset "NMRIO: pseudo-2D (Bruker pdata)" begin
    dat = exampledata("pseudo2D_XSTE")
    datC = loadnmr(artifact"pseudo2D_XSTE"; allcomponents=true)

    @test size(dat) == (2048, 10)

    @test parent(dat) isa Matrix{Float64}
    @test parent(datC) isa Matrix{Multicomplex{Float64,1,2}}

    @test dat[2, 3] == -32.08953857421875
    @test datC[2, 3] == Multicomplex(-32.08953857421875, -87.19561767578125)
    @test metadata(datC, 1, :mcindex) == 1 && isnothing(metadata(datC, 2, :mcindex))

    @test dat[:title] == "15N aSyn, 283 K - XSTE\nDELTA = 100ms, delta = 4ms, G = 5-95%"

    @test acqus(dat, :pulprog) == "stebpgp1s19xn.jk"

    @test dims(dat, 1) isa F1Dim
    @test dims(dat, 2) isa X2Dim
    @test dat[1, :window] == Cos2Window(0.0512)
    @test dat[1, :pseudodim] == false
    @test dat[2, :pseudodim] == true
end

@testset "pseudo3D data (Bruker pdata)" begin
    dat = exampledata("pseudo3D_HN_R2")
    @test size(dat) == (512, 512, 11)
    @test dat[5, 6, 7] == 2370.09375
    @test acqus(dat, :vclist) == [0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16]
end

@testset "3D data (Bruker pdata)" begin
    spec = exampledata("3D_HNCA")
    @test size(spec) == (432, 64, 256)
    @test spec[5, 6, 7] == 3506.28125
    @test label(spec, 1) == "1H"
    @test label(spec, 2) == "15N"
    @test label(spec, 3) == "13C"
    @test dims(spec, 1)[2] == 9.693887210893848
    @test dims(spec, 2)[end] == 107.84188874936109
    @test dims(spec, 3)[Near(45)] == 44.99702749174231

    spec13 = loadnmr(joinpath(artifact"3D_HNCA", "pdata", "131"))
    spec23 = loadnmr(joinpath(artifact"3D_HNCA", "pdata", "231"))
    @test size(spec13) == (512, 256)
    @test size(spec23) == (512, 128)
    @test spec23[Near(7.2), Near(122)] == 19735.930053710938
    @test spec13[Near(8.3), Near(63)] == -321.002197265625
    @test maximum(spec) / spec[:noise] ≈ 698.6861209993058
end

@testset "3D data (nmrPipe)" begin
    spec = loadnmr(joinpath(artifact"3D_HNCA_pipe", "smile", "test%03d.ft3"))
    @test size(spec) == (384, 512, 128)
    @test label(spec) == "500 uM 13C,15N ubiquitin, 298 K, BEST-HNCA"
    @test dims(spec, 1)[1] == 9.703234081887576
    @test dims(spec, 2)[1] == 68.7403353437212
    @test dims(spec, 3)[end] == 107.75223308061521
    @test spec[10, 11, 12] == 99473.6875
    @test spec[1, :window] == Cos2Window(0.09261055758378933)
    @test spec[2, :window] == CosWindow(0.057541115456864404)
    @test spec[3, :window] == CosWindow(0.03894846805778692)
    @test maximum(spec) / spec[:noise] == 492.00351859142154
end

@testset "NMRIO: 2D HN (ucsf)" begin
    dat = loadnmr(joinpath(artifact"2D_HN", "hmqc.ucsf"))

    @test size(dat) == (568, 256)
    @test length(dat) == 145408

    @test parent(dat) isa Matrix{Float64}

    @test dat[1, 100] == 59228.95703125
    @test dat[20, 100] == 27169.9140625
    @test dat[1] == -48482.26171875
    @test dat[1, end] == 36429.125
    @test dat[end, end] == -23958.6171875
    @test dat[end] == -23958.6171875
    @test maximum(dat) == 7.984743e6

    @test label(dat, 1) == "1H"
    @test label(dat, 2) == "15N"
end

@testset "3D data (ucsf)" begin
    spec = loadnmr(joinpath(artifact"3D_HNCA_ucsf", "hnca.ucsf"))
    @test size(spec) == (384, 512, 128)
    @test label(spec) == "500 uM 13C,15N ubiquitin, 298 K, BEST-HNCA"
    @test dims(spec, 1)[1] == 9.703234542832309
    @test dims(spec, 2)[1] == 68.74033212965945
    @test dims(spec, 3)[end] == 107.75222413893053
    @test spec[10, 11, 12] == 99473.6875
    @test spec[1, :window] isa UnknownWindow
    @test maximum(spec) / spec[:noise] ≈ 492.00351859142154
end

@testset "Samples and Experiments" begin
    spec = loadnmr(joinpath("test-data", "samples", "2"))
    @test spec[:date] == DateTime(2023, 10, 9, 10, 11, 34)
    @test basename(samplefile(spec)) == "2025-08-23_163924_lysozyme.json"
    @test sample(spec, "sample", "label") == "lysozyme"
    @test hassample(spec) == true
    @test sample(spec, :sample, :label) == "lysozyme"
    @test sample(spec, :sample, :components, :name) == ["HEWL"]

    spec = loadnmr(joinpath("test-data", "samples", "12"))
    @test sample(spec, "sample", "label") == "lysozyme (1mM Gd)"
    @test sample(spec, :sample, :components, :name) == ["HEWL", "gadodiamide"]

    spec = loadnmr(joinpath("test-data", "samples", "22"))
    @test isnothing(spec[:sample])
    @test hassample(spec) == false

    # Custom field auto-resolution
    s = NMRSample(joinpath("test-data", "samples", "custom_fields_fixture.json"))
    @test sample(s, :buffer, :solvent) == "95% H2O / 5% D2O + glycerol"
    @test sample(s, :sample, :components, :isotopic_labelling) ==
          ["2H,13C-Met", "unlabelled"]

    # NMRExperiment: same sample access as NMRData
    expt = NMRExperiment(joinpath("test-data", "samples", "2"))
    @test hassample(expt)
    @test sample(expt, :sample, :label) == "lysozyme"
    @test sample(expt, :sample, :components, :name) == ["HEWL"]

    expt2 = NMRExperiment(joinpath("test-data", "samples", "12"))
    @test sample(expt2, :sample, :components, :name) == ["HEWL", "gadodiamide"]

    # scansamples: returns sorted NMRSample list (excludes non-JSON files and experiment folders)
    samples = scansamples(joinpath("test-data", "samples"))
    @test length(samples) == 3   # lysozyme, lysozyme_1mM, custom_fields_fixture
    @test sample(samples[1], :sample, :label) == "lysozyme"
    @test sample(samples[2], :sample, :label) == "lysozyme (1mM Gd)"

    # findsample: match by timestamp from pre-scanned list
    @test sample(findsample(expt, samples), :sample, :label) == "lysozyme"
    @test sample(findsample(expt2, samples), :sample, :label) == "lysozyme (1mM Gd)"

    expt22 = NMRExperiment(joinpath("test-data", "samples", "22"))
    @test isnothing(findsample(expt22, samples))

    # scanexperiments: returns sorted NMRExperiment list, sample auto-associated
    expts = scanexperiments(joinpath("test-data", "samples"))
    @test length(expts) == 3
    @test expts[1][:date] == DateTime(2023, 10, 9, 10, 11, 34)
    @test expts[2][:date] == DateTime(2023, 10, 10, 11, 40, 17)
    @test expts[3][:date] == DateTime(2023, 10, 11, 10, 6, 7)
    @test label(expts[1]) == "10 mg/mL HEWL, 308 K, 1D zgesgp"
    @test label(expts[2]) == "10 mg/mL HEWL, 1 mM gadodiamide, 308 K, 1D zgesgp"
    @test label(expts[3]) == "10 mg/mL HEWL, 2 mM gadodiamide, 308 K, 1D zgesgp"
    @test hassample(expts[1])
    @test hassample(expts[2])
    @test !hassample(expts[3])
    @test sample(expts[1], :sample, :label) == "lysozyme"
    @test sample(expts[2], :sample, :label) == "lysozyme (1mM Gd)"

    # findsample from directory, and no-match case
    @test sample(findsample(expts[1], joinpath("test-data", "samples")), :sample, :label) ==
          "lysozyme"
    @test isnothing(findsample(expts[3]))

    # findexperiments: find all experiments belonging to a sample
    found = findexperiments(samples[1], expts)
    @test length(found) == 1
    @test found[1][:date] == DateTime(2023, 10, 9, 10, 11, 34)
end
@testset "NMRIO: peak lists" begin
    spec = exampledata("2D_HN")
    pks = detectpeaks(spec; sign=:positive)[1:5]
    pks[2] = Peak(pks[2]; label="G10")

    mktempdir() do dir
        # NMRAnalysis.jl peaklist.csv round trip
        f = joinpath(dir, "out", "peaklist.csv")
        @test savepeaks(f, pks) == f
        lines = readlines(f)
        @test startswith(lines[1], "#")
        @test lines[2] == "label,plane,x (ppm),y (ppm),lwx (Hz),lwy (Hz),intensity"
        @test startswith(lines[4], "G10,,")
        @test loadpeaks(f) == pks
        @test loadpeaks(f, spec) == pks
        @test_throws DimensionMismatch loadpeaks(f, exampledata("1D_19F"))

        # unlabelled peaks are written as X1, X2, ...
        savepeaks(f, [Peak(8.0, 120.0), Peak(7.0, 110.0; label="A")])
        @test label.(loadpeaks(f)) == ["X1", "A"]
        @test_throws ArgumentError savepeaks(f, [Peak(8.0, 120.0; label="A,B")])
        @test_throws ArgumentError savepeaks(f, [Peak(8.0, 1.0; label="A"),
                                                 Peak(7.0, 1.0; label="A")])

        # a 1D list
        f1 = joinpath(dir, "peaks1d.csv")
        p1 = detectpeaks(exampledata("1D_19F"))
        savepeaks(f1, p1)
        @test loadpeaks(f1) == p1

        # a hand-made list without a header
        fh = joinpath(dir, "hand.txt")
        write(fh, "# my peaks\nA 8.1 121.0\nB 7.5 118.2\n")
        hand = loadpeaks(fh)
        @test hand == [Peak(8.1, 121.0; label="A"), Peak(7.5, 118.2; label="B")]

        # a list written by NMRAnalysis.jl, with tracked positions and radius columns
        fn = joinpath(dir, "nmranalysis.csv")
        write(fn, """
              # NMRAnalysis.jl 0.4.3
              # X radius / ppm: 0.03
              label,plane,x (ppm),y (ppm),xradius (ppm),yradius (ppm)
              L23N,,8.5,122.1,0.03,0.3
              G10,1,8.1,110.2,0.03,0.3
              G10,2,8.12,110.4,0.03,0.3
              """)
        nmra = @test_logs (:warn,) loadpeaks(fn)
        @test label.(nmra) == ["L23N", "G10"]
        @test ppm(nmra[2]) == (8.1, 110.2)
        @test isnan(intensity(nmra[1]))

        # Sparky: written with the indirect dimension first, read back in NMRTools order
        fs = joinpath(dir, "hsqc.list")
        savepeaks(fs, pks)
        sl = readlines(fs)
        @test split(sl[1]) == ["Assignment", "w1", "w2", "Data", "Height"]
        @test parse(Float64, split(sl[4])[2]) ≈ ppm(pks[2], 2) atol = 1e-3
        sp = loadpeaks(fs)
        @test label.(sp) == label.(pks)
        @test all(isapprox.(ppm.(sp, 1), ppm.(pks, 1); atol=1e-3))
        @test all(isapprox.(ppm.(sp, 2), ppm.(pks, 2); atol=1e-3))

        # a Sparky list written the other way round is detected against the spectrum
        fsw = joinpath(dir, "swapped.list")
        write(fsw, "Assignment w1 w2\n\n" *
                   join(["?-? $(ppm(p, 1)) $(ppm(p, 2))" for p in pks], "\n"))
        sw = @test_logs (:info,) loadpeaks(fsw, spec)
        @test ppm.(sw, 1) ≈ ppm.(pks, 1)
        @test all(isempty, label.(sw))
    end
end

@testset "NMRIO: assigned shifts (NMR-STAR)" begin
    shifts = loadbmrb(joinpath("test-data", "bmrb", "bmr15000_shifts.str"))
    @test shifts isa Vector{ChemicalShift}
    @test length(shifts) == 340
    s = first(shifts)
    @test (s.residue, s.residuename, s.atom, nucleus(s), ppm(s)) == (2, "SER", "H", H1, 9.307)
    @test s.error == 0.01
    @test label(s) == "S2H"

    # backbone amides
    hn = peaklist(shifts, "H", "N")
    @test hn isa Vector{Peak{2}}
    @test length(hn) == 33
    @test hn[1] == Peak(9.307, 121.58; label="S2H-N")
    @test peaklist(shifts, "HN", "N") == hn
    @test all(p -> isnan(intensity(p)), hn)

    # inter-residue correlations
    con = peaklist(shifts, "CO-1", "N")
    @test con[1] == Peak(179.756, 121.409; label="D5C-F6N")
    hnco = peaklist(shifts, "H", "N", "C-1")
    @test hnco isa Vector{Peak{3}}
    @test hnco[1] == Peak(8.716, 121.409, 179.756; label="F6H-N-D5C")

    # all NH groups, including side chains
    nh = peaklist(shifts, H1, N15)
    @test length(nh) == 40
    @test label.(hn) ⊆ label.(nh)
    @test "W23HE1-NE1" in label.(nh)
    @test "R7HE-NE" in label.(nh)
    @test peaklist(shifts, N15, H1)[1] == Peak(121.58, 9.307; label="S2N-H")

    # CH groups: methylenes give one peak per proton, degenerate methyls just one
    ch = peaklist(shifts, H1, C13)
    l = label.(ch)
    @test "S2HB2-CB" in l && "S2HB3-CB" in l
    @test "V9HG2-CG2" in l
    @test !any(startswith("V9HG21"), l)
    @test all(p -> 0 < ppm(p, 1) < 12 && 5 < ppm(p, 2) < 180, ch)

    # a 1D list of every shift of one nucleus
    c = peaklist(shifts, C13)
    @test c isa Vector{Peak{1}}
    @test length(c) == count(s -> nucleus(s) == C13, shifts)
    @test "S2CB" in label.(c)

    # residue selection by filtering the shift list
    ser = filter(s -> s.residuename == "SER", shifts)
    @test all(p -> startswith(label(p), "S"), peaklist(ser, H1, C13))

    @test_throws ArgumentError peaklist(shifts, C13, N15)
    @test_throws ArgumentError peaklist(shifts, "C A")
    @test_throws ArgumentError loadbmrb(joinpath("test-data", "bmrb", "bmr15000_shifts.str");
                                          list=2)

    # shifts straight to a peak list
    f = joinpath("test-data", "bmrb", "bmr15000_shifts.str")
    @test loadbmrb(f, H1, N15) == nh
    @test loadbmrb(f, "H", "N", "CO-1") == hnco
    @test loadbmrb(f, C13) == c

    # the peak list can be saved for NMRAnalysis
    mktempdir() do dir
        f = joinpath(dir, "hsqc.csv")
        savepeaks(f, nh)
        @test loadpeaks(f) == nh
    end

    # NMR-STAR tokens: quoting, comments and multi-line text fields
    tokens = NMRTools.NMRIO._startokens("a 'b c' \"d'e\" it's # note\n;\nline 1\nline 2\n;\nf")
    @test tokens == ["a", "b c", "d'e", "it's", "\nline 1\nline 2", "f"]
end
