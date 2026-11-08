using NMRTools
using Documenter

DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRTools); recursive=true)
# fetch the example data first, so that download messages don't appear in doctest output
example("2D_HN")
doctest(NMRTools; manual=false)
