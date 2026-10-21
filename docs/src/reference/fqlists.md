# Frequency lists

Frequency lists are used in NMR experiments to specify arrays of frequencies for selective excitation, CEST, or other variable-frequency experiments. On Bruker spectrometers, these are stored in files like `fq1list` through `fq8list`.

NMRTools represents frequency lists as [`FQList`](@ref) structures that preserve encoding information (Hz vs ppm, relative vs absolute) and provide safe conversion functions.

## Encoding

Frequency lists can be specified on the spectrometer in several ways: in Hz or ppm, and relative to
the carrier or to the base frequency (0 ppm). An `FQList` stores the values with their unit and
whether they are relative, and, once attached to a channel with
[`withreference`](@ref NMRTools.NMRIO.withreference), the base and carrier frequencies needed to
convert them. The raw values are returned by `data(fqlist)`, but `ppm` and `hz` should be used to
read frequencies, since they apply the encoding.

## Accessing frequency lists

Frequency lists are automatically imported when loading data:

```julia
# Load CEST experiment
spec = loadnmr("path/to/cest/experiment")

# Get frequency list
fqlist = acqus(spec, :fq1list)
```

## Extracting frequency values

Use `ppm()` and `hz()` to safely extract frequency values. These functions handle all unit conversions automatically.

### Get chemical shifts (in ppm)

```julia
dim = dims(spec, F2Dim)
ppm_values = ppm(fqlist, dim)
```

### Get offsets from carrier (in Hz)

```julia
dim = dims(spec, F2Dim)
offset_hz = hz(fqlist, dim)
```

!!! tip
    Always use `ppm()` and `hz()` rather than accessing raw values directly. These functions handle all encoding schemes correctly.

## Creating frequency lists

Construct `FQList` objects programmatically:

```@example fqlist
using NMRTools

# Frequency list in Hz, relative to carrier
fq_rel = FQList([100.0, 200.0, 300.0], :Hz, true)

# Frequency list in ppm, absolute (chemical shifts)
fq_abs = FQList([7.5, 8.0, 8.5], :ppm, false)
```
