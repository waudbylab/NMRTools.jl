# Getting started

## Installing NMRTools

Install NMRTools.jl through the Julia package manager:

```julia
using Pkg
Pkg.add("NMRTools")
```

The examples in this tutorial also use the `CairoMakie` package for plotting, which can be obtained similarly.


## Plot a 1D spectrum

Let's load some example data. This can be a Bruker experiment directory, a specific pdata folder, or an NMRPipe-format file.

```@example 1
using NMRTools, CairoMakie
spec = exampledata("1D_19F")
```

NMRTools provides `nmrplot` for interactive and publication-quality figures. To plot a 1D spectrum:

```@example 1
fig = nmrplot(spec)
save("plot-y.svg", fig); nothing # hide
```

![](plot-y.svg)

We can zoom in using `xlims`, or select a chemical shift range from the data directly using the `..` selector — the y-axis then rescales to fit only the displayed region:

```@example 1
fig = nmrplot(spec[-124.5 .. -123])
save("plot-y2.svg", fig); nothing # hide
```

![](plot-y2.svg)

All plots can be saved as high quality vector graphics or PNG files using `save`:

```julia
save("myspectrum.pdf", fig)
save("myspectrum.svg", fig)
save("myspectrum.png", fig)
```


## Plot a 2D spectrum

Two-dimensional spectra are plotted in exactly the same way.

```@example 2d
using NMRTools, CairoMakie
spec = exampledata("2D_HN")
fig = nmrplot(spec)
save("plot-2d.svg", fig); nothing # hide
```

![](plot-2d.svg)

Contour levels start at five times the noise level with a geometric spacing of 1.7. Adjust the number of levels or spacing with the `levels` and `spacing` keywords. Use `title=""` to suppress the title:

```@example 2d
fig = nmrplot(spec; levels=4, spacing=2.0, title="")
save("plot-2d-scaled.svg", fig); nothing # hide
```

![](plot-2d-scaled.svg)


## Accessing your data

Spectrum data and associated axis information, metadata, etc, is encapsulated in an [`NMRData`](@ref) structure.
```@repl data
using NMRTools # hide
spec = exampledata("1D_19F")
```

Data can be accessed with conventional array indexing, but also using the value-based selectors, `Near` and `..`:
```@repl data
spec[100:105]
spec[Near(-124)]
spec[-124 .. -123.5]
```

This also works for multidimensional data. For example:
```@repl data
spec2d = exampledata("2D_HN")
spec2d[8.1 .. 8.3, Near(124)]
```

A plain array of data for the spectrum can be obtained from this using the [`data`](@ref) command:
```@repl data
data(spec)
```

Similarly, a plain vector containing axis values can be obtained from this using the [`data`](@ref) command, passing an additional argument to specify the dimension. This can either be a number or the axis type, e.g. `F1Dim`:
```@repl data
data(spec, 1)
```

## Next steps

- [Working with NMR data](@ref) covers loading, indexing, referencing, metadata and saving in more depth.
- [Processing raw data](@ref) shows how to turn raw `fid` and `ser` files into spectra, in one step with `autoprocess` or stage by stage:

```julia
spec = loadfid("myexperiment/1") |> autoprocess
```

- [Peaks and peak lists](@ref) covers peak picking and peak list files.
