# Peaks and peak lists

NMRTools represents a peak with the [`Peak`](@ref) type: a position in ppm, a height, a linewidth
in Hz and a label. A peak list is a `Vector{Peak{N}}`, so the usual tools for vectors (sorting,
filtering, indexing, `push!` and `deleteat!`) all work. This tutorial picks the peaks of a 2D
spectrum, looks at them, and saves them.

## Detecting peaks

[`detectpeaks`](@ref) finds the peaks in a 1D, 2D or higher-dimensional spectrum:

```@example peaks
using NMRTools, CairoMakie
spec = exampledata("2D_HN")
pks = detectpeaks(spec; sign=:positive)
```

A peak is a local maximum whose height exceeds a threshold, by default ten times the noise
level. Positions and heights are refined beyond the digital resolution by parabolic
interpolation, and the linewidth in each dimension is estimated from the half-height crossings.
Peaks are returned strongest first and labelled `X1`, `X2`, … to match the labels NMRAnalysis.jl
gives new peaks.

## Plotting peaks

`nmrplot!` overlays a peak list on a spectrum. 2D peaks are drawn as labelled crosses:

```@example peaks
fig, ax, _ = nmrplot(spec[7.5 .. 8.7, 116 .. 126])
nmrplot!(filter(p -> 7.5 < ppm(p, 1) < 8.7 && 116 < ppm(p, 2) < 126, pks))
save("peaks-2d.svg", fig); nothing # hide
```

![](peaks-2d.svg)

1D peaks are drawn as vertical lines, with their labels along the top of the axis:

```@example peaks
spec1d = exampledata("1D_19F")
fig, ax, _ = nmrplot(spec1d; xlims=(-123, -124.5))
nmrplot!(detectpeaks(spec1d))
save("peaks-1d.svg", fig); nothing # hide
```

![](peaks-1d.svg)

## Looking at a peak

Use the accessor functions to read a peak's properties:

```@example peaks
p = pks[1]
ppm(p)            # position, ordered like the spectrum's dimensions
```

```@example peaks
ppm(p, 2)         # one coordinate: here the ¹⁵N shift
```

```@example peaks
linewidth(p)      # full width at half maximum in each dimension, in Hz
```

`intensity(p)` gives the height and `label(p)` the label. Linewidths are kept in Hz because
they mean the same thing on every nucleus and at every field; the corresponding relaxation
rate is `π * linewidth(p, n)`. A quantity that is not known, such as the linewidth of a peak
read from a Sparky list, is `NaN`.

## Peak lists are tables

A peak list is a Tables.jl table, so it can be turned into a `DataFrame` or written with CSV.jl.
Sorting and filtering use the accessors:

```@example peaks
strong = sort(pks; by=intensity, rev=true)[1:5]
```

`ppm(n)` gives a function that reads coordinate `n` of a peak, so peaks can be sorted by position:

```@example peaks
lefttoright = sort(pks; by=ppm(1), rev=true)
ppm(1).(lefttoright[1:5])
```

```@example peaks
amides = filter(p -> 7.5 < ppm(p, 1) < 9, pks)
length(amides)
```

## Labelling peaks

Peaks are immutable, so changing a peak means making a new one. `Peak(p; ...)` copies a peak
with any of its fields replaced:

```@example peaks
pks[1] = Peak(pks[1]; label="G10")
pks[1]
```

[`findpeak`](@ref) finds the maximum nearest a given position, and labels it:

```@example peaks
findpeak(spec, (8.5, 120.0); label="L23")
```

## Saving the list

[`savepeaks`](@ref) writes a `peaklist.csv` that NMRAnalysis.jl reads directly, and
[`loadpeaks`](@ref) reads it back:

```@example peaks
file = joinpath(mktempdir(), "peaklist.csv")
savepeaks(file, pks)
loadpeaks(file)[1:3]
```

## Next steps

- [Assigning ubiquitin from the BMRB](bmrb.md): predict this spectrum's peaks from assigned
  chemical shifts, and register them onto it
- [Detect peaks](../howto/detectpeaks.md): thresholds, negative peaks, regions and 3D data
- [Save and load peak lists](../howto/peaklists.md): file formats, including Sparky
