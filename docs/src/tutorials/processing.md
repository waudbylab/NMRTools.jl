# Processing a 1D spectrum

!!! warning "Under development"
    Loading and processing raw time-domain data is new and its interface may change.

This tutorial turns a raw ¹H FID into a spectrum, first in one step with [`autoprocess`](@ref) and
then one step at a time, plotting the data after each step. The FID is a 1D ¹H spectrum of
ubiquitin in H₂O, recorded with excitation sculpting to suppress the water signal.

```@example proc
using NMRTools, CairoMakie
fid = exampledata("1D_1H"; fid=true)
nothing # hide
```

For your own data, use [`loadfid`](@ref NMRTools.NMRIO.loadfid) with the experiment folder:
`fid = loadfid("path/to/experiment/1")`.

## In one step

```@example proc
spec = autoprocess(fid)
fig = nmrplot(spec[-1 .. 10])
save("proc-auto.svg", fig); nothing # hide
```

![](proc-auto.svg)

`autoprocess` removes the residual water signal with a solvent filter (for ¹H spectra in H₂O),
applies a cos² window, zero-fills to twice the size, Fourier transforms and phases the spectrum
automatically, and returns the real part. The rest of this page goes through these steps,
leaving the water in; see [Solvent filtering](@ref) for the filter.

## The FID

`loadfid` returns time-domain data, with a time axis in seconds:

```@example proc
dims(fid)
```

```@example proc
fig = nmrplot(fid)
save("proc-fid.svg", fig); nothing # hide
```

![](proc-fid.svg)

The FID is complex. `nmrplot` draws the real part; the imaginary part is `imag.(fid)`. Zooming in
on the start shows the two parts, and the points before t = 0:

```@example proc
start = fid[-0.007 .. 0.03]
fig, ax, _ = nmrplot(start; title="start of the FID")
nmrplot!(ax, imag.(start))
vlines!(ax, [0.0]; color=:grey, linestyle=:dash)
save("proc-fid-start.svg", fig); nothing # hide
```

![](proc-fid-start.svg)

The points before t = 0 are the build-up of the spectrometer's digital filter. `loadfid` keeps
them, and the Fourier transform accounts for them; see
[Bruker data and TopSpin conventions](../concepts/bruker.md#Digital-filtering-and-group-delay).

## Apodization

[`apodize`](@ref) multiplies the FID by a window function, here a cos² window that falls to zero
at the end of the acquisition time. This removes the truncation artefacts ("sinc wiggles") that
an FID cut off before it has decayed would give:

```@example proc
windowed = apodize(fid, Cos2Window())
fig, ax, _ = nmrplot(fid; title="before and after apodization")
nmrplot!(ax, windowed)
save("proc-apodize.svg", fig); nothing # hide
```

![](proc-apodize.svg)

The window is recorded in the metadata of the time dimension:

```@example proc
windowed[T1Dim, :window]
```

See [Window functions](../reference/windows.md) for the other windows.

## Zero-filling

[`zerofill`](@ref) appends zeros to the FID, doubling its length by default. This gives more
points in the spectrum, which interpolates the lineshapes:

```@example proc
filled = zerofill(windowed)
size(windowed), size(filled)
```

```@example proc
fig = nmrplot(filled; title="zero-filled")
save("proc-zerofill.svg", fig); nothing # hide
```

![](proc-zerofill.svg)

## Fourier transform

[`ft`](@ref) turns the time dimension into a frequency dimension:

```@example proc
unphased = ft(filled)
dims(unphased)
```

The spectrum is still complex, and not yet phased:

```@example proc
fig = nmrplot(real(unphased)[-1 .. 10]; title="after Fourier transform")
save("proc-ft.svg", fig); nothing # hide
```

![](proc-ft.svg)

## Phasing

[`autophase`](@ref) estimates the zero- and first-order phase corrections and applies them:

```@example proc
phased = autophase(unphased)
phased[F1Dim, :phase]
```

The corrections can also be set by hand with [`phase`](@ref), in degrees. The values use the same
convention as TopSpin's `PHC0` and `PHC1`, so phases from TopSpin can be copied directly:

```@example proc
manual = phase(unphased, F1Dim, 16.8, 0.0)
fig, ax, _ = nmrplot(real(phased)[6 .. 10]; title="amide region")
nmrplot!(ax, real(manual)[6 .. 10])
save("proc-phase.svg", fig); nothing # hide
```

![](proc-phase.svg)

## The real spectrum

The phased spectrum still has an imaginary part. `real` takes the real part, giving a
`Float64` spectrum:

```@example proc
final = real(phased)
eltype(phased), eltype(final)
```

```@example proc
fig = nmrplot(final[-1 .. 10])
save("proc-final.svg", fig); nothing # hide
```

![](proc-final.svg)

## As a pipeline

Each processing function returns a new spectrum, and called without data it returns a function,
so the steps can be chained with `|>`:

```@example proc
spec = fid |> apodize(Cos2Window()) |> zerofill |> ft |> autophase |> real
nothing # hide
```

## Next steps

- [Process step by step](../howto/process.md): window functions, [solvent filtering](@ref "Solvent filtering"),
  phasing options, extracting regions and pseudo-2D data
- [Process a 2D with linear prediction](../howto/linearprediction.md)
- [Save spectra and FIDs](../howto/save.md)
- [Bruker data and TopSpin conventions](../concepts/bruker.md)
