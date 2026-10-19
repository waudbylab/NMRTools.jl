# Peaks and peak lists

NMRTools represents a peak with the `Peak{N}` type: a position in ppm, a height, a linewidth
in Hz and a label. A peak list is simply a `Vector{Peak{N}}`, so the usual tools for vectors
(sorting, filtering, indexing, `push!` and `deleteat!`) all work, and a peak list is also a
Tables.jl table, so `DataFrame(pks)` works directly.

Peak lists are written in the `peaklist.csv` format used by
[NMRAnalysis.jl](https://github.com/waudbylab/NMRAnalysis.jl), so peaks picked here can be
loaded straight into an NMRAnalysis analysis, and vice versa.

## Detecting peaks

`detectpeaks` finds the peaks in a 1D, 2D or higher-dimensional spectrum:

```@example peaks
using NMRTools, CairoMakie
spec = exampledata("2D_HN")
pks = detectpeaks(spec; sign=:positive)
```

A peak is a local maximum whose height exceeds a threshold, by default ten times the noise
level. Positions and heights are refined beyond the digital resolution by parabolic
interpolation, and the linewidth in each dimension is estimated from the half-height
crossings. Peaks are returned strongest first and labelled `X1`, `X2`, … to match the
labels NMRAnalysis gives new peaks.

The threshold can be set relative to the noise or absolutely, and negative peaks can be
included or picked on their own:

```julia
detectpeaks(spec; snr=20)              # only peaks above 20σ
detectpeaks(spec; threshold=1e6)       # an absolute threshold
detectpeaks(spec; sign=:negative)      # negative peaks only; the default is :both
```

To pick peaks in a pseudo-2D experiment, index one plane first, e.g.
`detectpeaks(spec[:, 1])`.

Like the processing functions, `detectpeaks`, `register` and `savepeaks` have forms for use in a
pipeline:

```julia
pks = loadfid("hsqc/1") |> autoprocess |> realest |> detectpeaks(snr=20)
pks |> savepeaks("peaklist.csv")
loadbmrb(15000, H1, N15) |> register(spec) |> savepeaks("bmrb.csv")
```

## Plotting peaks

`nmrplot!` overlays a peak list on a spectrum. 2D peaks are drawn as labelled crosses:

```@example peaks
fig, ax, _ = nmrplot(spec)
nmrplot!(pks)
fig
```

and 1D peaks as vertical lines, with their labels along the top of the axis:

```@example peaks
spec1d = exampledata("1D_19F")
fig, ax, _ = nmrplot(spec1d; xlims=(-123, -124.5))
nmrplot!(detectpeaks(spec1d))
fig
```

With Plots.jl, `plot!(pks)` does the same.

## Working with peaks

Use the accessor functions rather than reaching into the fields:

```@example peaks
p = pks[1]
ppm(p)            # position, ordered like the spectrum's dimensions
```

```@example peaks
linewidth(p)      # full width at half maximum in each dimension, in Hz
```

`intensity(p)` gives the height and `label(p)` the label. Linewidths are kept in Hz because
they mean the same thing on every nucleus and at every field; the corresponding relaxation
rate is `π * linewidth(p, n)`. A quantity that is not known, such as the linewidth of a peak
read from a Sparky list, is `NaN`.

Peaks are immutable, so changing a peak means making a new one. `Peak(p; ...)` copies a peak
with any of its fields replaced, which makes assignment straightforward:

```@example peaks
pks[1] = Peak(pks[1]; label="G10")
pks[1]
```

Peaks can also be created by hand, or snapped onto the nearest maximum in a spectrum with
`findpeak`:

```@example peaks
findpeak(spec, (8.5, 120.0); label="L23")
```

Sorting and filtering use the same accessors:

```julia
sort(pks; by=intensity, rev=true)           # strongest first
sort(pks; by=p -> ppm(p, 1), rev=true)      # left to right
filter(p -> 7.5 < ppm(p, 1) < 9, pks)       # a region of the spectrum
```

## Saving and loading

`savepeaks` writes a `peaklist.csv` that NMRAnalysis.jl reads directly, and `loadpeaks`
reads it back, along with lists written by NMRAnalysis and simple hand-made lists of
`label x y` lines:

```julia
savepeaks("peaklist.csv", pks)
pks = loadpeaks("peaklist.csv")
```

The file has one row per peak, with units in the column names:

```
# NMRTools.jl peak list
label,plane,x (ppm),y (ppm),lwx (Hz),lwy (Hz),intensity
G10,,8.5013,120.371,23.1,14.2,2.53e7
```

`x` is the direct dimension and `y` the indirect one. The blank `plane` column means the
position applies to every plane of a series, which is how NMRAnalysis reads it. Labels
identify peaks in NMRAnalysis, so they must be unique; unlabelled peaks are written as `X1`,
`X2`, ….

Sparky `.list` files are written and read too. Sparky lists the indirect dimension first,
which is the reverse of NMRTools' order. Passing the spectrum to `loadpeaks` checks the
order against where the shifts actually fall, so a list written the other way round is
still read correctly:

```julia
savepeaks("hsqc.list", pks)
pks = loadpeaks("hsqc.list", spec)
```

## Registering a peak list onto a spectrum

A peak list made for one spectrum rarely lines up exactly with another, because of differences in
referencing, temperature or sample conditions. This applies to a list carried over from an earlier
spectrum of the same sample, and even more to one predicted from assigned chemical shifts (see
[Assigned chemical shifts from the BMRB](@ref)). `register!` finds the single offset that best
matches the list to the observed peaks and moves the list by it:

```julia
register!(pks, spec)                        # reports the offset applied
register!(pks, spec; maxshift=(0.3, 2.0))   # limit the search, in ppm per dimension
```

Every offset within `maxshift` (by default 10% of the spectral width) is tried, scoring how
much signal lies under the shifted peaks. Absolute intensities are used, so negative peaks
count as much as positive ones, and the score grows only logarithmically with intensity, so
a single very strong peak cannot outweigh many weaker matches. `register` returns a moved
copy instead. Only a uniform translation is applied; individual peaks are not moved.
