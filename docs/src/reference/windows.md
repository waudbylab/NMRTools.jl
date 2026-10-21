# Window functions

Window functions (also called apodization functions) are applied to time-domain NMR data before Fourier transformation to modify the lineshape and signal-to-noise characteristics of the resulting spectrum. NMRTools automatically detects and stores the window function used during processing, whether by TopSpin or NMRPipe for loaded spectra, or by [`apodize`](@ref) for raw data processed in NMRTools (see [Processing a 1D spectrum](../tutorials/processing.md)). When a window is passed to `apodize` without `tmax`, for example `Cos2Window()`, its length is set from the acquisition time of each dimension.

## Overview

Window functions are represented by subtypes of the abstract type [`WindowFunction`](@ref). Each window function type stores:
- Parameters specific to that window function
- The acquisition time `tmax` (calculated after linear prediction but before zero filling)

Window functions are stored as axis metadata and can be accessed through the `:window` key:

```@example windows
using NMRTools

# Load example spectrum
spec = exampledata("2D_HN")

# Get the window function for a specific dimension
window = spec[F2Dim, :window]
```

## Available window functions

`t` is the time from t = 0 and `tmax` the acquisition time. Omit `tmax` when passing a window to
[`apodize`](@ref) and it is set from each dimension's acquisition time.

| Window | w(t) | Notes |
|:-------|:-----|:------|
| `NullWindow(tmax)` | ``1`` | no apodization |
| `UnknownWindow(tmax)` | – | the window could not be determined from the processing parameters |
| `ExponentialWindow(lb, tmax)` | ``\exp(-\pi \, \mathrm{lb} \, t)`` | line broadening `lb` (Hz): improves S/N at the cost of resolution; negative `lb` narrows lines |
| `CosWindow(tmax)` | ``\cos\left(\frac{\pi t}{2 t_\mathrm{max}}\right)`` | cosine (sine bell shifted by 90°); common in indirect dimensions |
| `Cos2Window(tmax)` | ``\cos^2\left(\frac{\pi t}{2 t_\mathrm{max}}\right)`` | squared cosine; falls to zero more smoothly than `CosWindow`. Alias `Cos²Window` |
| `GeneralSineWindow(offset, endpoint, power, tmax)` | ``\sin^{\mathrm{power}}\left(\pi\,\mathrm{offset} + \frac{(\mathrm{endpoint} - \mathrm{offset})\,\pi t}{t_\mathrm{max}}\right)`` | general sine bell; `offset` and `endpoint` (0 to 1) set the phase at the start and end. `offset=0.5, endpoint=1` gives `CosWindow` (`power=1`) or `Cos2Window` (`power=2`) |
| `LorentzToGaussWindow(expHz, gaussHz, tmax)` | ``\exp(\pi \, \mathrm{expHz} \, t) \exp\left(-\frac{(\pi \, \mathrm{gaussHz} \, t)^2}{2}\right)`` | Lorentz-to-Gauss transformation: removes `expHz` (Hz) of Lorentzian broadening and adds `gaussHz` (Hz) of Gaussian broadening; resolution enhancement |
| `GeneralGaussWindow(expHz, gaussHz, center, tmax)` | ``\exp(\pi \, \mathrm{expHz} \, t) \exp\left(-\frac{(\pi \, \mathrm{gaussHz} \, (t - t_\mathrm{center}))^2}{2}\right)`` | as `LorentzToGaussWindow`, with the Gaussian maximum at `center` (0 to 1) of the acquisition time; `center=0` gives `LorentzToGaussWindow` |

```@example windows
ExponentialWindow(5.0, 0.1), Cos2Window(0.1), GeneralSineWindow(0.5, 1.0, 2.0, 0.1)
```

## Working with window functions

### Extracting window parameters

Window function objects store their parameters as fields:

```julia
# Load spectrum
spec = exampledata("2D_HN")

# Get window for F2 dimension
win = spec[F2Dim, :window]

# For ExponentialWindow
if win isa ExponentialWindow
    println("Line broadening: ", win.lb, " Hz")
    println("Acquisition time: ", win.tmax, " s")
end

# For sine windows
if win isa CosWindow
    println("Acquisition time: ", win.tmax, " s")
end
```

### Computing the apodization function

The [`apod`](@ref) function returns the time-domain apodization vector:

```@example windows
# Load a spectrum
spec = exampledata("2D_HN")

# Get apodization function for F2 dimension
apod_function = apod(spec, F2Dim)

println("Apodization vector length: ", length(apod_function))
println("First few values: ", apod_function[1:5])
```

You can also compute it directly from a dimension:

```julia
dim = dims(spec, F2Dim)
apod_function = apod(dim)
```

By default, the apodization vector is zero-filled to match `tdzf`. To get only the acquired points:

```julia
apod_function = apod(spec, F2Dim, false)  # zerofill=false
```

### Simulating lineshapes

The [`lineshape`](@ref) function generates a simulated spectrum including the effect of the window function:

```@example windows
# Get the ¹H axis
axis = dims(spec, F1Dim)

# Simulate a peak at 8.5 ppm with R2 = 10 s-1
chemical_shift = 8.5  # ppm
R2 = 10.0            # s-1
simulated = lineshape(axis, chemical_shift, R2)

println("Simulated lineshape length: ", length(simulated))
```

This is useful for:
- Fitting experimental peaks
- Understanding window function effects
- Quality control and validation

You can also request complex-valued lineshapes:

```julia
# Get complex lineshape
simulated_complex = lineshape(axis, chemical_shift, R2, ComplexLineshape())
```
