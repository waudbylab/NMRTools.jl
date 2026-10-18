using Documenter, NMRTools, Plots

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRTools); recursive=true)

makedocs(;
         modules=[NMRTools],
         format=Documenter.HTML(; prettyurls=(get(ENV, "CI", nothing) == "true"),
                                canonical="https://waudbylab.org/NMRTools.jl/stable/",
                                # all API reference is on one page (api.md)
                                size_threshold_warn=300 * 2^10,
                                size_threshold=500 * 2^10),
         pages=["Home" => "index.md",
                "Getting started" => "quickstart.md",
                "Tutorials" => ["Working with NMR data" => "tutorials/data.md",
                                "Processing raw data" => "tutorials/processing.md",
                                "Plotting with Plots" => "tutorials/plotting.md",
                                "Plotting with Makie" => "tutorials/makie.md",
                                "Peaks and peak lists" => "tutorials/peaks.md"],
                "Worked examples" => ["1D relaxation analysis" => "tutorials/relaxation.md",
                                      "Diffusion analysis" => "tutorials/diffusion.md",
                                      "Animating spectra" => "tutorials/animation.md"],
                "Metadata" => ["Overview" => "metadata/index.md",
                               "Times and frequencies" => "metadata/timeunits.md",
                               "Power levels" => "metadata/power.md",
                               "Frequency lists" => "metadata/fqlists.md",
                               "Window functions" => "metadata/windows.md",
                               "Sample metadata" => "metadata/samples.md"],
                "Utilities" => "utilities.md",
                "Reference guide" => ["Chemical shift referencing" => "reference/referencing.md",
                                      "Metadata" => "reference/metadata.md",
                                      "Annotations" => "reference/annotations.md",
                                      "Coherences and isotope data" => "reference/coherences.md"],
                "API" => "api.md",
                "Index" => "indexes.md"],
         # repo=GitHub("waudbylab", "NMRTools.jl"),
         sitename="NMRTools.jl",
         authors="Chris Waudby",
         warnonly=[:missing_docs],)

deploydocs(;
           repo="github.com/waudbylab/NMRTools.jl.git",
           devbranch="main",)
