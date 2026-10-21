using NMRTools
using Documenter

DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRTools); recursive=true)
doctest(NMRTools; manual=false)
