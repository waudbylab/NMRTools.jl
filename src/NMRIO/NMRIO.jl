module NMRIO

using ..NMRBase
import ..NMRBase: NMRSample, NMRExperiment
using Artifacts
using ArtifactUtils
using Dates
using LazyArtifacts
using JSON
using MulticomplexNumbers
using YAML

include("../sample-schemas/migration-code/migrate.jl")
using .SchemaMigrate

include("lists.jl")
include("loadnmr.jl")
include("jdx.jl")
include("acqus.jl")
include("nmrpipe.jl")
include("bruker.jl")
include("annotation.jl")
include("samples.jl")
include("experiment.jl")
include("scan.jl")
include("ucsf.jl")
include("exampledata.jl")
include("sumexpts.jl")

export loadnmr
export loadmetadata
export loadjdx
export FQList
export sumexpts
export scanexperiments
export scansamples
export findsample
export findexperiments

# examples
export exampledata

end
