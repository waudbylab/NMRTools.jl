using NMRTools
using Test
using SafeTestsets

@safetestset "NMRBase" begin
    include("nmrbase_test.jl")
end

@safetestset "NMRIO" begin
    include("nmrio_test.jl")
end

@safetestset "Processing" begin
    include("processing_test.jl")
end

@safetestset "Advanced processing" begin
    include("advanced_processing_test.jl")
end

@safetestset "loadfid" begin
    include("loadfid_test.jl")
end

@safetestset "Autophase" begin
    include("autophase_test.jl")
end

@safetestset "Bruker writing" begin
    include("brukerwrite_test.jl")
end

@safetestset "PlotsExt" begin
    include("plotsext_test.jl")
end

@safetestset "MakieExt" begin
    include("makieext_test.jl")
end

@safetestset "Doctests" begin
    include("doctest_test.jl")
end
