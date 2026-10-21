# Makie plot options

Keyword arguments of [`nmrplot`](@ref) for Makie. Makie's own `Axis` attributes can be passed through the `axis` keyword.

## 1D spectra

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `title` | spectrum label | Pass `""` to suppress. |
| `color` | auto (Wong palette) | Line colour. |
| `normalize` | `true` | Scale by scans and receiver gain. Pass `false` to use raw intensities. |
| `vstack` | `false` | Stack spectra vertically. Pass a number to scale the offset. |
| `xlims` | auto | Chemical shift range (ppm). y-axis rescales to fit. |
| `legend` | `false` | Legend position, e.g. `:rt`. Entries named from spectrum labels. |
| `axis` | `(;)` | Extra `Axis` keyword arguments (escape hatch). |

## 2D spectra

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `title` | spectrum label | Pass `""` to suppress. |
| `color` / `poscolor` | auto (Wong palette) | Positive contour colour. |
| `negcontours` | `true` | Show negative contours. |
| `negcolor` | faded positive | Negative contour colour. |
| `normalize` | `true` | Scale contour threshold by scans/gain. Pass a reference spectrum for cross-spectrum comparison. |
| `levels` | `12` | Number of contour levels (Int), or a vector of explicit positive level values. |
| `spacing` | `1.7` | Geometric ratio between successive contour levels; also the step size per ↑/↓ keypress. |
| `xprojection` | `nothing` | 1D projection along the direct dimension. Pass `true` for max-projection or a 1D `NMRData`. |
| `yprojection` | `nothing` | 1D projection along the indirect dimension. |
| `xlims` | auto | Direct-dimension range (ppm). |
| `ylims` | auto | Indirect-dimension range (ppm). |
| `legend` | `false` | Legend position, e.g. `:lt`. Entries named from spectrum labels. |
| `axis` | `(;)` | Extra `Axis` keyword arguments (escape hatch). |

## Pseudo-2D spectra

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `style` | `:heatmap` | `:heatmap`, `:flat` (overlaid 1D slices) or `:waterfall` (3D slice view) |
