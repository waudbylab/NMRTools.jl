# Metadata overview

[`NMRData`](@ref) objects contain metadata on processing and acquisition parameters, populated automatically when loading a spectrum.

## Types of metadata

Metadata is organised into two categories:

1. **Spectrum metadata**: information about the experiment as a whole: number of scans, receiver gain, pulse programme, title, temperature, the noise level (calculated when the spectrum is loaded), acquisition parameters from the `acqus` file, and the contents of auxiliary files such as vclists and vdlists.
2. **Axis metadata**: information specific to each dimension: number of points, the original time-domain size (before zero filling, linear prediction or extraction of a subregion), carrier frequency, spectral width, window function and phase correction.

Metadata are stored in dictionaries keyed by symbols such as `:ns` or `:pulseprogram`. The standard keys are listed in the [Metadata](../reference/metadata.md) reference page.

## Quick start

### Accessing spectrum metadata

```@example overview
using NMRTools

# Load example spectrum
spec = exampledata("2D_HN")

# The whole dictionary
metadata(spec)
```

```@example overview
# A single entry, using the metadata function
metadata(spec, :ns)
```

```@example overview
# Or use dictionary-style lookup
spec[:pulseprogram]
```

### Accessing axis metadata

```@example overview
# By dimension type
spec[F2Dim, :label]
```

```@example overview
# Or by dimension number
spec[2, :swhz]
```

```@example overview
# All metadata for one axis
metadata(spec, F1Dim)
```

### Accessing acquisition parameters

```@example overview
# Get acquisition parameter
acqus(spec, :te)
```

```@example overview
# Access arrayed parameters
acqus(spec, :p, 1)  # First pulse length
```

Parameter names can be given as strings (case insensitive, e.g. `"TE"`) or as lowercase symbols (e.g. `:te`).

!!! note
    Arrayed parameters such as pulse lengths and delays are returned as dictionaries rather than lists, to avoid confusion between Bruker arrays, which are zero-based, and Julia arrays, which are one-based.

### Spectrometer channels

Frequency dimensions carry their own `:bf`, `:sf` and `:nucleus` metadata, but these only describe nuclei that were *detected*, i.e. that have an axis. To reach a nucleus that was manipulated but never detected, for example the carbon channel in a ¹³C-edited experiment read out on ¹H, use the spectrometer **channel** model.

Channels describe the physical RF channels (`:f1`…`:f8`) independently of the detected dimensions. `channels(spec)` returns a dictionary of all channels, and `channels(spec, ...)` returns a single channel by label, by `Nucleus`, or by a string giving either:

```@example overview
channels(spec, :f1)        # by spectrometer channel label
```

```@example overview
channels(spec, "13C")[:bf] # carbon base frequency, even with no carbon axis
```

Each channel dictionary has keys `:label`, `:nucleus`, `:bf` (base frequency, Hz), `:sfo` (carrier frequency, Hz), `:offsethz` and `:offsetppm`.

## Special metadata types

NMRTools provides specialized types for handling complex metadata:

### Window functions

Window functions (apodization) applied during processing are stored as structured objects that preserve all parameters:

```@example overview
# Get window function for a dimension
win = spec[F2Dim, :window]
```

See the [Window functions](windows.md) page for detailed documentation.

### Power levels

Power levels are represented as [`Power`](@ref) objects that handle conversion between Watts and dB:

```julia
# Power levels from acqus file
p = acqus(spec, :pl, 1)
db(p)      # Get as dB
watts(p)   # Get as Watts
```

See the [Power levels](power.md) page for detailed documentation.

### Frequency lists

Frequency lists (from `fq1list`, etc.) are stored as [`FQList`](@ref) objects that preserve unit and reference information:

```julia
# Get frequency list
fqlist = acqus(spec, :fq1list)

# Extract as ppm or Hz
ppm(fqlist, dims(spec, F2Dim))
hz(fqlist, dims(spec, F2Dim))
```

See the [Frequency lists](fqlists.md) page for detailed documentation.

### Sample metadata

Sample metadata describing the physical sample, buffer composition, and experimental details can be automatically loaded when present:

```julia
# Check if sample metadata is available
hassample(spec)

# Access sample metadata
sample(spec)

# Navigate nested metadata
sample(spec, :sample, :label)
```

See the [Sample metadata](samples.md) page for detailed documentation.

## Auxiliary files

Files such as `vclist`, `vdlist`, and frequency lists are automatically imported when present:

```@example overview
# Access variable delay list
relaxation_spec = exampledata("pseudo3D_HN_R2")
acqus(relaxation_spec, :vclist)
```

NMRTools performs automatic unit conversion:
- `:vclist`: Variable loop counter (dimensionless)
- `:vdlist`: Variable delays, in seconds
- `:valist`: Variable amplitude, in dB (converted from Watts if necessary)
- `:vplist`: Variable pulse lengths, in seconds
- `:fq1list` through `:fq8list`: Frequency lists (see [Frequency lists](fqlists.md))

## Next steps

- **[Times and frequencies](timeunits.md)**: Units used for pulse lengths, delays and frequencies
- **[Power levels](power.md)**: Representation of powers and RF field strength calculations
- **[Frequency lists](fqlists.md)**: Frequency list handling and conversions
- **[Window functions](windows.md)**: Apodization functions and lineshape effects
- **[Sample metadata](samples.md)**: Working with sample metadata and the NMRSamples schema
- **[Reference: Metadata](../reference/metadata.md)**: metadata functions and the standard metadata keys
