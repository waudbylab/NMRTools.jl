# Detect peaks

[`detectpeaks`](@ref) finds the peaks in a 1D, 2D or higher-dimensional spectrum and returns them
as a `Vector{Peak{N}}`, strongest first. For an introduction, see the tutorial
[Peaks and peak lists](../tutorials/peaks.md).

```@example detect
using NMRTools, CairoMakie
spec = exampledata("2D_HN")
nothing # hide
```

## Thresholds

A peak is a local maximum whose height exceeds a threshold, by default ten times the noise level
(`spec[:noise]`). Set the threshold relative to the noise with `snr`, or absolutely with
`threshold`:

```@example detect
length(detectpeaks(spec)), length(detectpeaks(spec; snr=30))
```

```julia
detectpeaks(spec; threshold=1e6)      # an absolute threshold
```

## Positive and negative peaks

By default both positive and negative peaks are found. Select one sign with `sign`:

```@example detect
pos = detectpeaks(spec; sign=:positive)
neg = detectpeaks(spec; sign=:negative)
length(pos), length(neg)
```

```@example detect
fig, ax, _ = nmrplot(spec; title="positive (black) and negative (red) peaks")
nmrplot!(ax, pos)
nmrplot!(ax, neg; color=:red)
save("detect-sign.svg", fig); nothing # hide
```

![](detect-sign.svg)

## A region of the spectrum

Select a region before picking:

```@example detect
amides = detectpeaks(spec[7.5 .. 8.7, 116 .. 126]; sign=:positive)
fig, ax, _ = nmrplot(spec[7.5 .. 8.7, 116 .. 126])
nmrplot!(ax, amides)
save("detect-region.svg", fig); nothing # hide
```

![](detect-region.svg)

or filter the list afterwards by position:

```julia
filter(p -> 7.5 < ppm(p, 1) < 8.7, pks)
```

## 1D and 3D spectra

The same function works in any number of dimensions:

```@example detect
detectpeaks(exampledata("1D_19F"))
```

```@example detect
hnca = detectpeaks(exampledata("3D_HNCA"); snr=20, sign=:positive)
hnca[1:3]
```

For a pseudo-2D series, pick peaks in one plane: `detectpeaks(spec[:, 1])`.

## A single peak

[`findpeak`](@ref) finds the maximum nearest a given position, refined by interpolation, and labels
it:

```@example detect
findpeak(spec, (8.27, 121.0); label="A28")
```

[`findnearest`](@ref) gives the index of the peak in a list nearest a position. `scale` weights the
dimensions, since a ppm in ¹⁵N is not a ppm in ¹H:

```@example detect
i = findnearest(pos, (8.27, 121.0); scale=(1, 0.2))
pos[i]
```

## Filtering by properties

A peak list is a vector, so `filter` and `sort` with the accessors [`intensity`](@ref),
[`linewidth`](@ref), `ppm` and `label` select peaks:

```@example detect
strong = filter(p -> intensity(p) > 50 * spec[:noise], pos)
narrow = filter(p -> linewidth(p, 1) < 30, pos)
length(strong), length(narrow)
```

## In a pipeline

```julia
pks = loadfid("hsqc/1") |> autoprocess |> realest |> detectpeaks(snr=20)
```
