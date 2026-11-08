# Getting started

This page takes about ten minutes. It loads and plots a 1D and a 2D spectrum, and processes a raw
FID into a spectrum.

## Installing NMRTools

Install NMRTools.jl and the CairoMakie plotting package through the Julia package manager:

```julia
using Pkg
Pkg.add(["NMRTools", "CairoMakie"])
```

Use GLMakie instead of CairoMakie for interactive plots in a window.

## Plot a 1D spectrum

[`loadnmr`](@ref NMRTools.NMRIO.loadnmr) loads a spectrum from the path of a Bruker experiment
folder or an nmrPipe file. Example data, from
[NMRExampleData.jl](https://waudbylab.org/NMRExampleData.jl/stable/), are found with
`example(name)`, which gives the path to load, downloading the dataset the first time.

```@example gs
using NMRTools, CairoMakie
spec = loadnmr(example("1D_19F"))
size(spec), label(spec)
```

[`nmrplot`](@ref) plots it, with axes labelled from the spectrum's metadata:

```@example gs
fig = nmrplot(spec)
save("gs-1d.svg", fig); nothing # hide
```

![](gs-1d.svg)

Select a chemical shift range with `..` to zoom in. The y-axis rescales to fit:

```@example gs
fig = nmrplot(spec[-124.5 .. -123])
save("gs-1d-zoom.svg", fig); nothing # hide
```

![](gs-1d-zoom.svg)

## Plot a 2D spectrum

Two-dimensional spectra are plotted in the same way. Contour levels start at five times the noise
level:

```@example gs
spec2d = loadnmr(example("2d_hn_ubq"))
fig = nmrplot(spec2d)
save("gs-2d.svg", fig); nothing # hide
```

![](gs-2d.svg)

Ranges can be selected in both dimensions, and `levels` and `spacing` change the contours:

```@example gs
fig = nmrplot(spec2d[7 .. 9, 115 .. 125]; levels=8, spacing=1.5, title="")
save("gs-2d-zoom.svg", fig); nothing # hide
```

![](gs-2d-zoom.svg)

Figures can be saved as vector graphics or bitmaps with `save("spectrum.pdf", fig)`.

## Look inside a spectrum

A spectrum is an array with axes and metadata attached:

```@example gs
spec2d[Near(8.5), Near(120)]       # the intensity nearest 8.5 ppm, 120 ppm
```

```@example gs
spec2d[:ns], spec2d[F2Dim, :nucleus]   # number of scans; nucleus of the second dimension
```

## Process a FID

[`loadfid`](@ref NMRTools.NMRIO.loadfid) loads raw time-domain data, and
[`autoprocess`](@ref) turns it into a phased spectrum: apodization, zero-filling, Fourier
transformation and automatic phasing.

```@example gs
fid = loadfid(example("1d_1h_ubq"))
spec1h = autoprocess(fid)
fig = nmrplot(spec1h[-1 .. 10]; title="")
save("gs-fid.svg", fig); nothing # hide
```

![](gs-fid.svg)

`autoprocess` returns the real spectrum. For a ¹H spectrum in H₂O, it also removes the residual
water signal with a solvent filter.

## Next steps

- [Working with spectra](tutorials/spectra.md): indexing, axes, metadata and saving
- [Processing a 1D spectrum](tutorials/processing.md): what `autoprocess` does, step by step
- [Peaks and peak lists](tutorials/peaks.md)
- [Assigning ubiquitin from the BMRB](tutorials/bmrb.md)
- Tools: [combining experiments](tools/sumexpts.md) and
  [peak lists from the BMRB](tools/bmrb.md)
