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

## Peak lists from BMRB assignments

`loadshifts` reads the assigned chemical shifts of a BMRB entry, downloading it by its ID,
or from a local NMR-STAR file. `peaklist` then assembles them into the peaks expected in a
given experiment:

```julia
shifts = loadshifts(15000)               # download BMRB entry 15000

hsqc = peaklist(shifts, "H", "N")        # backbone amide ¹H,¹⁵N HSQC
nmrplot(hsqc; xlabel="¹H (ppm)", ylabel="¹⁵N (ppm)")
```

Dimensions are given in spectrum order, direct dimension first. Each argument names an atom,
optionally with an offset to a neighbouring residue, so any triple-resonance correlation can
be written down directly:

```julia
peaklist(shifts, "CO-1", "N")            # CON: C′ of residue i-1 with N of residue i
peaklist(shifts, "H", "N", "CO-1")       # HNCO
peaklist(shifts, "H", "N", "CA")         # HNCA, i peaks
peaklist(shifts, "H", "N", "CA-1")       # HNCA, i-1 peaks
```

Passing nuclei instead of atom names gives every directly bonded X–H pair, which is what an
HSQC actually shows. This includes side-chain NH and NH₂ groups for ¹⁵N, and every CH, CH₂
and CH₃ group for ¹³C. Methyl and methylene protons with identical shifts give a single peak:

```julia
peaklist(shifts, H1, N15)                # ¹H,¹⁵N HSQC including side chains
peaklist(shifts, H1, C13)                # ¹H,¹³C HSQC
```

A shift list is a `Vector{ChemicalShift}`, so residues can be selected with `filter` before
building peaks, e.g. for the methyl region of an ILV-labelled sample:

```julia
ilv = filter(s -> s.residuename in ("ILE", "LEU", "VAL"), shifts)
peaklist(ilv, H1, C13)
```

Peaks are labelled by residue and atom (`L23H-N`, `D5C-F6N`, `V9HG2-CG2`), so a predicted
list can be overlaid on a spectrum with `nmrplot!`, or saved with `savepeaks` as the starting
point for an analysis in NMRAnalysis.jl. Remember that BMRB shifts were measured under the
deposited sample conditions and referencing, so expect systematic offsets.
