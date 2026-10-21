# Working with spectra

NMR measurements are arrays of data with coordinates attached to each axis: chemical shifts,
times, gradient strengths. NMRTools stores them as [`NMRData`](@ref) structures, which hold the
data, the axes and metadata on acquisition and processing. This tutorial loads some spectra,
takes them apart, and plots each piece along the way.

```@example spectra
using NMRTools, CairoMakie
```

## Loading a spectrum

[`loadnmr`](@ref NMRTools.NMRIO.loadnmr) loads processed data: a Bruker experiment folder, a
`pdata` folder, an nmrPipe file or a UCSF file. Here we use the bundled example data:

```@example spectra
spec = exampledata("1D_19F")
nothing # hide
```

`exampledata(name; path=true)` gives the folder the data were loaded from, which is what you would
pass to `loadnmr` for your own data:

```julia
spec = loadnmr("path/to/experiment/1")             # pdata/1 by default
spec = loadnmr("path/to/experiment/1/pdata/101")   # another processed spectrum
```

See [Load spectra](../howto/load.md) for the other formats.

## Indexing and selecting

An `NMRData` can be indexed like any Julia array:

```@example spectra
spec[1000:1004]
```

It can also be indexed by chemical shift, with the selectors `Near(x)` (the nearest point),
`At(x)` (an exact match) and `a .. b` (everything between `a` and `b`):

```@example spectra
spec[Near(-123.7)]
```

Selecting a range returns another `NMRData`, with its axis cut to match, so it can be plotted
straight away:

```@example spectra
region = spec[-124.5 .. -123]
fig = nmrplot(region)
save("spectra-region.svg", fig); nothing # hide
```

![](spectra-region.svg)

The same works in two dimensions. A range in both dimensions gives a smaller 2D spectrum:

```@example spectra
spec2d = exampledata("2D_HN")
fig = nmrplot(spec2d[7.5 .. 8.7, 116 .. 126])
save("spectra-2d-region.svg", fig); nothing # hide
```

![](spectra-2d-region.svg)

while a single point in one dimension gives a 1D slice through the other:

```@example spectra
slice = spec2d[:, Near(121.0)]
fig = nmrplot(slice; title="¹H trace at 121 ppm ¹⁵N")
save("spectra-2d-slice.svg", fig); nothing # hide
```

![](spectra-2d-slice.svg)

!!! warning "Shared metadata"
    Slices share the metadata dictionary of the spectrum they came from, so changing the metadata
    of one changes the other. Use `deepcopy` to make an independent copy.

## Axes

The axes of a spectrum are its dimensions:

```@example spectra
dims(spec2d)
```

Each has a type that says what it measures: here `F1Dim` and `F2Dim`, the direct and indirect
frequency dimensions. A dimension can be used like a vector of its coordinates, and
[`data`](@ref) gives the plain numbers:

```@example spectra
data(spec2d, F2Dim)[1:5]
```

`data(spec2d)` returns the plain matrix of intensities. See [Dimensions](../concepts/dimensions.md)
for the full set of dimension types.

## Series of spectra

A pseudo-2D experiment, such as a relaxation series, has a frequency dimension and a
non-frequency dimension:

```@example spectra
relax = exampledata("pseudo2D_T2")
dims(relax)
```

The second dimension is loaded as `X2Dim`, an unknown non-frequency dimension. Once we know the
relaxation delays (here 4 ms per loop of the vclist), [`setrelaxtimes`](@ref) replaces it with a
`TrelaxDim`, and plots label it correctly:

```@example spectra
relax = setrelaxtimes(relax, acqus(relax, :vclist) * 4e-3, "s")
fig = nmrplot(relax[0.7 .. 1.0, :]; style=:waterfall)
save("spectra-pseudo2d.svg", fig); nothing # hide
```

![](spectra-pseudo2d.svg)

Each column is a 1D spectrum:

```@example spectra
fig = nmrplot(relax[0.7 .. 1.0, 1])
nmrplot!(fig, relax[0.7 .. 1.0, end])
save("spectra-pseudo2d-slices.svg", fig); nothing # hide
```

![](spectra-pseudo2d-slices.svg)

## Real and imaginary parts

Processed spectra are loaded as real data by default. Pass `allcomponents=true` to load the
imaginary parts too:

```@example spectra
spec19f = loadnmr(exampledata("1D_19F"; path=true); allcomponents=true)
eltype(spec19f)
```

The element type is a first-order multicomplex number: a complex number whose imaginary unit is
`im1`. `real` and `imag` give the two parts as real spectra:

```@example spectra
region = spec19f[-124.5 .. -123]
fig, ax, _ = nmrplot(real.(region); title="real and imaginary parts")
nmrplot!(ax, imag.(region))
save("spectra-realimag.svg", fig); nothing # hide
```

![](spectra-realimag.svg)

Multidimensional data have an imaginary unit for each dimension; see
[Multicomplex numbers and quadrature](../concepts/multicomplex.md).

## Arithmetic

Arithmetic works on `NMRData` as on any array, and keeps the axes and metadata:

```@example spectra
halved = spec2d - 0.5 * spec2d
fig = nmrplot(halved[7.5 .. 8.7, 116 .. 126]; title="spec2d - 0.5 * spec2d")
save("spectra-arith.svg", fig); nothing # hide
```

![](spectra-arith.svg)

## Metadata

Metadata describe the experiment as a whole and each of its axes. Look up spectrum metadata with
a key, and axis metadata with a dimension and a key:

```@example spectra
spec2d[:ns], spec2d[:pulseprogram], spec2d[:temperature]
```

```@example spectra
spec2d[F1Dim, :nucleus], spec2d[F1Dim, :swppm]
```

Acquisition parameters from the `acqus` files are read with [`acqus`](@ref):

```@example spectra
acqus(spec2d, :te), acqus(spec2d, :p, 1)
```

Times are in seconds and frequencies in Hz. See [Metadata](../concepts/metadata.md) for how
metadata are organised, and [Find acquisition parameters](../howto/acqus.md) for more.

## Saving

A processed spectrum can be written back as Bruker processed data, to view in TopSpin:

```julia
savenmr("path/to/experiment/1", spec; procno=101)    # writes pdata/101
```

See [Save spectra and FIDs](../howto/save.md).

## Next steps

- [Processing a 1D spectrum](processing.md): from a raw FID to a spectrum
- [Plot 1D spectra](../howto/plot1d.md) and [Plot 2D spectra](../howto/plot2d.md)
- [The NMRData model](../concepts/nmrdata.md)
