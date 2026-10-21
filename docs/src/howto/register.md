# Register a peak list onto a spectrum

A peak list made for one spectrum rarely lines up exactly with another, because of differences in
referencing, temperature or sample conditions. This applies to a list carried over from an earlier
spectrum of the same sample, and even more to one predicted from assigned chemical shifts.
[`register!`](@ref) finds the single offset that best matches the list to the observed peaks, and
moves the list by it. `register` returns a moved copy instead.

```@example reg
using NMRTools, CairoMakie
spec = exampledata("2D_HN")
nothing # hide
```

## A list from the BMRB

```@example reg
predicted = peaklist(exampledata("BMRB_ubiquitin"), "H", "N")
moved = register(predicted, spec)
nothing # hide
```

The offset applied is logged. See [Assigning ubiquitin from the BMRB](../tutorials/bmrb.md) for a
full example.

## A list from another spectrum

A list picked in one spectrum can be registered onto another, for example the next point of a
titration:

```@example reg
series = exampledata("2D_HN_titration")
pks = detectpeaks(series[1]; snr=20, sign=:positive)
moved = register(pks, series[end])
fig, ax, _ = nmrplot(series[end][7.5 .. 8.7, 116 .. 126]; title="")
nmrplot!(ax, filter(p -> 7.5 < ppm(p, 1) < 8.7 && 116 < ppm(p, 2) < 126, moved))
save("reg-series.svg", fig); nothing # hide
```

![](reg-series.svg)

## A saved list

Load the list with [`loadpeaks`](@ref) first:

```julia
pks = loadpeaks("peaklist.csv")
register!(pks, spec)
```

## Limiting the search

Every offset within `maxshift` (by default 10% of the spectral width) is tried. Give `maxshift` in
ppm, one value per dimension, to narrow the search, or `0` to hold a dimension fixed:

```julia
register!(pks, spec; maxshift=(0.3, 2.0))   # at most ±0.3 ppm ¹H, ±2 ppm ¹⁵N
register!(pks, spec; maxshift=(0.3, 0))     # shift ¹H only
```

## How the offset is chosen

Each candidate offset is scored by how much signal lies under the shifted peaks. Absolute
intensities are used, so negative peaks count as much as positive ones, and the score grows only
logarithmically with intensity, so a single strong peak cannot outweigh many weaker matches. The
best offset is refined to a fraction of a point.

Only a uniform translation is applied; individual peaks are not moved. Use [`findpeak`](@ref) on
each registered peak to snap it to the nearest maximum:

```@example reg
snapped = [findpeak(series[end], ppm(p); label=label(p)) for p in moved[1:3]]
```

## Higher dimensions

Registration works in any number of dimensions, as long as the list and spectrum have the same
number:

```julia
register!(peaklist(shifts, "H", "N", "CA"), exampledata("3D_HNCA"))
```

## In a pipeline

```julia
loadbmrb(17769, H1, N15) |> register(spec) |> savepeaks("bmrb.csv")
```
