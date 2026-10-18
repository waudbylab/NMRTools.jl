module NMRIO

using ..NMRBase
import ..NMRBase: NMRSample, NMRExperiment
using Artifacts
using Dates
using Downloads
using LazyArtifacts
using JSON
using MulticomplexNumbers
using YAML

include("../sample-schemas/migration-code/migrate.jl")
using .SchemaMigrate

include("lists.jl")
include("loadnmr.jl")
include("loadfid.jl")
include("brukerwrite.jl")
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
include("peaklists.jl")
include("bmrb.jl")

export loadnmr
export loadfid
export savefid, savenmr
export loadmetadata
export loadjdx
export FQList
export withreference
export sumexpts
export loadpeaks, savepeaks
export loadbmrb
export scanexperiments
export scansamples
export findsample
export findexperiments

# examples
export exampledata

end
