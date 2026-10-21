using Documenter, NMRTools, Plots

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRTools); recursive=true)

# Fetch the example data before building: on a fresh machine the first use of each artifact
# prints a download message, which would end up in example and doctest output.
foreach(name -> exampledata(name; path=true), NMRTools.NMRIO.EXAMPLEDATA_NAMES)

# the release notes page is generated from NEWS.md at the repository root
write(joinpath(@__DIR__, "src", "releasenotes.md"),
      replace(read(joinpath(@__DIR__, "..", "NEWS.md"), String),
              r"^# NMRTools.jl release notes" => "# Release notes"))

makedocs(;
         modules=[NMRTools],
         format=Documenter.HTML(; prettyurls=(get(ENV, "CI", nothing) == "true"),
                                canonical="https://waudbylab.org/NMRTools.jl/stable/"),
         pages=["Home" => "index.md",
                "Getting started" => "getting-started.md",
                "Tutorials" => ["Working with spectra" => "tutorials/spectra.md",
                                "Processing a 1D spectrum" => "tutorials/processing.md",
                                "Peaks and peak lists" => "tutorials/peaks.md",
                                "Assigning ubiquitin from the BMRB" => "tutorials/bmrb.md",
                                "Relaxation analysis" => "tutorials/relaxation.md"],
                "Tools" => ["Combining experiments" => "tools/sumexpts.md",
                            "Peak lists from the BMRB" => "tools/bmrb.md"],
                "How-to guides" => ["Data" => ["Load spectra" => "howto/load.md",
                                               "Process step by step" => "howto/process.md",
                                               "Process a 2D with linear prediction" => "howto/linearprediction.md",
                                               "Reference chemical shifts" => "howto/referencing.md",
                                               "Save spectra and FIDs" => "howto/save.md"],
                                    "Plotting" => ["Plot 1D spectra" => "howto/plot1d.md",
                                                   "Plot 2D spectra" => "howto/plot2d.md",
                                                   "Animate plots" => "howto/animate.md",
                                                   "Plot with Plots.jl" => "howto/plots.md"],
                                    "Peaks" => ["Detect peaks" => "howto/detectpeaks.md",
                                                "Register a peak list onto a spectrum" => "howto/register.md",
                                                "Save and load peak lists" => "howto/peaklists.md"],
                                    "Metadata" => ["Find acquisition parameters" => "howto/acqus.md",
                                                   "Find sample information" => "howto/samples.md"]],
                "Concepts" => ["The NMRData model" => "concepts/nmrdata.md",
                               "Dimensions" => "concepts/dimensions.md",
                               "Multicomplex numbers and quadrature" => "concepts/multicomplex.md",
                               "Metadata" => "concepts/metadata.md",
                               "Time and frequency units" => "concepts/units.md",
                               "Bruker data and TopSpin conventions" => "concepts/bruker.md"],
                "Reference" => ["Metadata keys" => "reference/metadatakeys.md",
                                "Window functions" => "reference/windows.md",
                                "Power levels" => "reference/power.md",
                                "Frequency lists" => "reference/fqlists.md",
                                "Sample metadata schema" => "reference/samples.md",
                                "Nuclei, coherences and Xi ratios" => "reference/nuclei.md",
                                "Makie plot options" => "reference/makieoptions.md"],
                "Experimental" => ["Pulse programme annotations" => "experimental/annotations.md"],
                "API" => ["Loading and saving" => "api/io.md",
                          "Processing" => "api/processing.md",
                          "Peaks and assignments" => "api/peaks.md",
                          "Plotting" => "api/plotting.md",
                          "Dimensions" => "api/dimensions.md",
                          "Metadata" => "api/metadata.md",
                          "Window functions and lineshapes" => "api/windows.md",
                          "Nuclei and coherences" => "api/nuclei.md",
                          "Internals" => "api/internals.md"],
                "Release notes" => "releasenotes.md"],
         sitename="NMRTools.jl",
         authors="Chris Waudby",
         checkdocs=:exports,
         warnonly=[:missing_docs],)

# `@repl` blocks show errors instead of failing the build: fail if any page outside the API
# reference (whose doctests show intended errors) contains one
let build = joinpath(@__DIR__, "build")
    failed = String[]
    for (root, _, files) in walkdir(build), f in files
        endswith(f, ".html") || continue
        path = joinpath(root, f)
        startswith(relpath(path, build), "api") && continue
        occursin("ERROR: ", read(path, String)) && push!(failed, relpath(path, build))
    end
    isempty(failed) || error("errors in the output of example blocks on: $(join(failed, ", "))")
end

deploydocs(;
           repo="github.com/waudbylab/NMRTools.jl.git",
           devbranch="main",)
