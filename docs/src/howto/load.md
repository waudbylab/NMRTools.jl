# Load spectra

[`loadnmr`](@ref NMRTools.NMRIO.loadnmr) loads processed spectra and detects the format from the
path. [`loadfid`](@ref NMRTools.NMRIO.loadfid) loads raw time-domain data.

```@example load
using NMRTools
```

## Bruker processed data

Pass an experiment folder to load `pdata/1`, or a `pdata` folder to load another processed
spectrum:

```julia
spec = loadnmr("path/to/experiment/1")               # pdata/1
spec = loadnmr("path/to/experiment/1/pdata/101")     # pdata/101
```

```@example load
dir = exampledata("2D_HN"; path=true)
spec = loadnmr(dir)
size(spec)
```

1D, 2D and 3D processed data (`1r`, `2rr`, `3rrr`) are read. By default only the real part is
loaded; `allcomponents=true` loads the imaginary parts as well, as multicomplex data:

```@example load
spec = loadnmr(dir; allcomponents=true)
eltype(spec)
```

## Raw data

`loadfid` reads a `fid` (1D) or `ser` (nD) file. Pass the experiment folder or the file itself:

```@example load
fid = loadfid(dir)
size(fid), eltype(fid)
```

On loading, indirect dimensions are recombined from their quadrature scheme, the digital filter's
group delay is recorded, the first point is scaled, and the data are scaled to match TopSpin's
intensities. See [Bruker data and TopSpin conventions](../concepts/bruker.md#What-loadfid-does) for
the details, and [Process step by step](process.md) for what to do next.

## nmrPipe

Pass the file name of an nmrPipe spectrum. For a 3D stored as a series of planes, give a template
with `%03d` (or similar) in place of the plane number:

```@example load
pipe = exampledata("3D_HNCA_pipe"; path=true)
basename(pipe)
```

```@example load
spec3d = loadnmr(pipe)
size(spec3d)
```

1D, 2D and 3D nmrPipe data are read, including pseudo-3D plane series.

## UCSF (Sparky)

```@example load
ucsf = loadnmr(exampledata("3D_HNCA_ucsf"; path=true))
size(ucsf)
```

## Attaching acquisition metadata

`loadnmr` looks for the Bruker acquisition files (`acqus`, pulse programme, lists) next to the
spectrum, so that acquisition metadata are available for nmrPipe or UCSF data processed in the
experiment folder. If the spectrum was saved elsewhere, give the experiment folder with
`experimentfolder`:

```julia
spec = loadnmr("processed/hsqc.ft2"; experimentfolder="data/hsqc/1")
```

## Series of experiments

Load each experiment with `loadnmr` and collect them in a vector. Plotting functions accept vectors
of spectra:

```julia
spectra = [loadnmr("titration/$i") for i in 1:11]
```

[`scanexperiments`](@ref NMRTools.NMRIO.scanexperiments) lists the experiments in a folder without
loading their data; see [Find sample information](samples.md).

## Noise level

When a spectrum is loaded its noise level is estimated and stored as `:noise`. Plots use it to set
contour levels, and [`detectpeaks`](@ref) to set thresholds:

```@example load
spec[:noise]
```

[`estimatenoise`](@ref) re-estimates it, for example after editing the data.
