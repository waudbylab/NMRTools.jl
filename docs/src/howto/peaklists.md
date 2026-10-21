# Save and load peak lists

[`savepeaks`](@ref) writes a peak list to a file and [`loadpeaks`](@ref) reads one. The format is
chosen from the file extension.

```@example pl
using NMRTools
spec = exampledata("2D_HN")
pks = detectpeaks(spec; snr=20, sign=:positive)
dir = mktempdir()
nothing # hide
```

## NMRAnalysis.jl format (`.csv`)

`savepeaks` writes a `peaklist.csv` that [NMRAnalysis.jl](https://waudbylab.org/NMRAnalysis.jl)
reads directly:

```@example pl
file = joinpath(dir, "peaklist.csv")
savepeaks(file, pks)
print(join(readlines(file)[1:4], "\n"))
```

The file has one row per peak, with units in the column names. `x` is the direct dimension and `y`
the indirect one. The blank `plane` column means the position applies to every plane of a series,
which is how NMRAnalysis reads it. Labels identify peaks in NMRAnalysis, so they must be unique;
unlabelled peaks are written as `X1`, `X2`, ….

`loadpeaks` reads it back, along with lists written by NMRAnalysis and simple hand-made lists of
`label x y` lines:

```@example pl
loadpeaks(file)[1:3]
```

## Sparky format (`.list`)

Sparky lists the indirect dimension first, the reverse of NMRTools' order. Passing the spectrum to
`loadpeaks` checks the order against where the shifts fall, so a list written the other way round
is still read correctly:

```@example pl
sparky = joinpath(dir, "hsqc.list")
savepeaks(sparky, pks)
print(join(readlines(sparky)[1:4], "\n"))
```

```@example pl
loadpeaks(sparky, spec)[1:3]
```

Quantities that Sparky lists don't record, such as linewidths, are `NaN` after loading.

## As tables

A peak list is a Tables.jl table, so it converts to a `DataFrame` or any other table type:

```julia
using DataFrames
df = DataFrame(pks)
```

## In a pipeline

```julia
pks |> savepeaks("peaklist.csv")
detectpeaks(spec) |> savepeaks("peaklist.csv")
```
