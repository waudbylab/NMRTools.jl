# Window functions

Window functions (also called apodization functions) are applied to time-domain NMR data before Fourier transformation to modify the lineshape and signal-to-noise characteristics of the resulting spectrum. NMRTools automatically detects and stores the window function used during processing, whether by TopSpin or NMRPipe for loaded spectra, or by [`apodize`](@ref) for raw data processed in NMRTools (see [Processing raw data](@ref)). When a window is passed to `apodize` without `tmax`, for example `Cos2Window()`, its length is set from the acquisition time of each dimension.

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

### No apodization

**`NullWindow(tmax)`**

Represents no apodization applied to the data. The multiplication factor is 1.0 for all time points.

```@example windows
w = NullWindow(0.1)  # 0.1 second acquisition time
```

### Unknown apodization

**`UnknownWindow(tmax)`**

Used when the window function cannot be determined from the processing parameters. This typically occurs with older data formats or non-standard processing.

### Exponential (line broadening)

**`ExponentialWindow(lb, tmax)`**

Exponential apodization with line broadening `lb` in Hz. This is the most common apodization function in NMR processing.

```math
w(t) = \exp(-\pi \cdot \mathrm{lb} \cdot t)
```

```@example windows
# 5 Hz line broadening, 0.1 s acquisition
w = ExponentialWindow(5.0, 0.1)
```

**Effect:**
- Positive `lb`: Increases line broadening, improves S/N, reduces resolution
- Negative `lb`: Decreases line broadening (resolution enhancement), degrades S/N
- Typical values: 1-10 Hz for routine 1D proton spectra

**Example:**
```@example windows
# Get line broadening from a loaded spectrum
lb = spec[F2Dim, :window]
println("Line broadening: ", lb)
```

### Sine/Cosine windows

**`CosWindow(tmax)`**

Pure cosine (half-period) window function:

```math
w(t) = \cos\left(\frac{\pi t}{2t_{\mathrm{max}}}\right)
```

```@example windows
w = CosWindow(0.1)
```

**Effect:** Commonly used for indirect dimensions in multidimensional NMR. Improves resolution with moderate S/N penalty.

---

**`Cos2Window(tmax)`**

Squared cosine window function:

```math
w(t) = \cos^2\left(\frac{\pi t}{2t_{\mathrm{max}}}\right)
```

```@example windows
w = Cos2Window(0.1)
```

**Effect:** Gentler apodization than cosine, better S/N but slightly lower resolution.

---

**`GeneralSineWindow(offset, endpoint, power, tmax)`**

General sine window function with adjustable parameters:

```math
w(t) = \left[\sin\left(
    \pi \cdot \mathrm{offset} +
    \frac{(\mathrm{endpoint} - \mathrm{offset})\pi t}{t_{\mathrm{max}}}
    \right)\right]^{\mathrm{power}}
```

**Parameters:**
- `offset`: Initial phase (0 to 1), determines starting value
- `endpoint`: Final phase (0 to 1), determines ending value
- `power`: Exponent applied to sine function
- `tmax`: Acquisition time

```@example windows
# Typical shifted sine-bell
w = GeneralSineWindow(0.5, 1.0, 2.0, 0.1)
```

**Special cases:**
- `offset=0.5, power=1.0` → `CosWindow`
- `offset=0.5, power=2.0` → `Cos2Window`

### Gaussian windows

**`LorentzToGaussWindow(expHz, gaussHz, tmax)`**

Lorentz-to-Gauss transformation, combining inverse exponential with Gaussian broadening:

```math
w(t) = \exp(+\pi \cdot \mathrm{expHz} \cdot t) \cdot \exp\left(-\frac{(\pi \cdot \mathrm{gaussHz} \cdot t)^2}{2}\right)
```

```@example windows
# 5 Hz inverse exponential, 10 Hz Gaussian
w = LorentzToGaussWindow(5.0, 10.0, 0.1)
```

**Parameters:**
- `expHz`: Inverse exponential (line narrowing) in Hz
- `gaussHz`: Gaussian broadening in Hz
- `tmax`: Acquisition time

**Effect:** Converts Lorentzian lineshapes (natural in NMR) to Gaussian lineshapes. Improves resolution for crowded spectra but requires careful parameter optimization.

---

**`GeneralGaussWindow(expHz, gaussHz, center, tmax)`**

General Gaussian window with adjustable center position:

```math
w(t) = \exp(+\pi \cdot \mathrm{expHz} \cdot t) \cdot
       \exp\left(-\frac{(\pi \cdot \mathrm{gaussHz} \cdot (t - t_{\mathrm{center}}))^2}{2}\right)
```

**Parameters:**
- `expHz`: Inverse exponential in Hz
- `gaussHz`: Gaussian broadening in Hz
- `center`: Position of Gaussian maximum (0 to 1)
- `tmax`: Acquisition time

```@example windows
# Shifted Gaussian window
w = GeneralGaussWindow(5.0, 10.0, 0.3, 0.1)
```

**When center = 0:** Simplifies to `LorentzToGaussWindow`

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
# Get the F2 axis
axis = dims(spec, F2Dim)

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
