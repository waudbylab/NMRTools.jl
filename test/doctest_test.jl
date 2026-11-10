using NMRTools
using Documenter

DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRTools); recursive=true)
# fetch the example data first, so that download messages don't appear in doctest output
example("2d_hn_ubq")
doctest(NMRTools; manual=false)
