using NMRTools
using Documenter

DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRTools); recursive=true)
# fetch the example data first, so that download messages don't appear in doctest output
exampledata("2D_HN"; path=true)
doctest(NMRTools; manual=false)
