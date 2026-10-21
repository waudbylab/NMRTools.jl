# Plot 2D spectra

[`nmrplot`](@ref) plots 2D spectra as contour plots, with contour levels set from the noise level.

```@example 1
using NMRTools
using CairoMakie   # static figures; use GLMakie for interactive
```

## A single spectrum

```@example 1
spec2d = exampledata("2D_HN")
fig = nmrplot(spec2d)
save("makie-2D.svg", fig); nothing # hide
```

![](makie-2D.svg)

Add 1D projections along either axis:

```@example 1
fig = nmrplot(exampledata("2D_HN"), xprojection=true, yprojection=true)
save("makie-2D-projections.svg", fig); nothing # hide
```

![](makie-2D-projections.svg)

Plot a series of spectra in one call; colours cycle automatically and a legend can be added:

```@example 1
fig = nmrplot(exampledata("2D_HN_titration"), legend=:lt)
save("makie-2D-series.svg", fig); nothing # hide
```

![](makie-2D-series.svg)

## Contour colours and levels

Negative contours are shown by default as a faded version of the positive colour. Disable them or set a different colour:

```@example 1
fig = nmrplot(spec2d; negcontours=false)
save("makie-2Db.svg", fig); nothing # hide
```

![](makie-2Db.svg)

```@example 1
fig = nmrplot(spec2d; negcolor=:red)
save("makie-2Dc.svg", fig); nothing # hide
```

![](makie-2Dc.svg)

Pass `levels` as an integer to change the number of contour levels, or as a vector of absolute values to use explicit levels:

```julia
nmrplot(spec2d; levels=8)                    # 8 geometric levels from 5σ
nmrplot(spec2d; levels=[2e6, 4e6, 8e6])      # explicit positive levels
```

`spacing` sets the geometric ratio between levels (default 1.7):

```julia
nmrplot(spec2d; spacing=2.0, levels=10)
```

## Keyboard navigation (GLMakie)

In GLMakie, press **↑** / **↓** to raise or lower all contour levels in the figure by one step (multiplied/divided by `spacing`). The key repeats while held. This works for any spectrum added with `nmrplot` or `nmrplot!`.

## Pseudo-2D spectra

Pseudo-2D spectra (one frequency + one non-frequency dimension) support three display styles via the `style` keyword:

| `style` | Description |
|:--------|:------------|
| `:heatmap` | colour map (default) |
| `:flat` | overlaid 1D slices |
| `:waterfall` | 3D slice view |

```@example 1
diffusiondata = exampledata("pseudo2D_XSTE")
diffusiondata = setgradientlist(diffusiondata, LinRange(0.02, 0.98, 10))
fig = nmrplot(diffusiondata, xlims=(7, 9.5))
save("makie-pseudo2D-heatmap.svg", fig); nothing # hide
```

![](makie-pseudo2D-heatmap.svg)

```@example 1
fig = nmrplot(diffusiondata, xlims=(7, 9.5), style=:flat)
save("makie-pseudo2D-flat.svg", fig); nothing # hide
```

![](makie-pseudo2D-flat.svg)

```@example 1
fig = nmrplot(diffusiondata, xlims=(7, 9.5), style=:waterfall)
save("makie-pseudo2D-waterfall.svg", fig); nothing # hide
```

![](makie-pseudo2D-waterfall.svg)


## Normalising a series

Pass a reference spectrum as `normalize` to put several spectra on the same contour levels, taking
account of differences in scans and receiver gain:

```@example 1
spectra = exampledata("2D_HN_titration")
fig = nmrplot(spectra[[1, 11]]; normalize=spectra[1], title="first and last titration points")
save("makie-2D-normalize.svg", fig); nothing # hide
```

![](makie-2D-normalize.svg)

## 3D spectra

Plot a plane or a strip of a 3D spectrum by selecting it first:

```@example 1
hnca = exampledata("3D_HNCA")
fig = nmrplot(hnca[:, Near(120.0), :]; title="¹H,¹³C plane at 120 ppm ¹⁵N")
save("makie-3D-plane.svg", fig); nothing # hide
```

![](makie-3D-plane.svg)

## Peaks

Overlay a peak list with `nmrplot!(ax, peaks)`; see [Peaks and peak lists](../tutorials/peaks.md).

## Saving figures and layout

Saving and composing figures works as for 1D spectra; see [Plot 1D spectra](plot1d.md#Saving-figures).

See [Makie plot options](../reference/makieoptions.md) for every keyword.
