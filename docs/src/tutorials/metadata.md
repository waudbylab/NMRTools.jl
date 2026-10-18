# Metadata

[`NMRData`](@ref) objects contain metadata on processing and acquisition parameters, populated automatically when loading a spectrum. Entries are divided into **spectrum metadata** – associated with the experiment in general – and **axis metadata**, that are associated with a particular dimension of the data.

Examples of spectrum metadata include: number of scans, receiver gain, pulse program, experiment title, number of dimensions, noise level (calculated when the spectrum is loaded), acquisition parameters (pulse lengths etc, from the `acqus` file), contents of auxiliary files (e.g. vclists and vdlists).

Examples of axis metadata include: the number of points, the original time domain size (before zero filling, linear prediction or extraction of a subregion), carrier frequency, spectrum width, window function used for processing.


## Accessing spectrum metadata

Metadata are stored in a dictionary labelled by symbols such as `:ns` or `:pulseprogram`. This dictionary can be accessed using the [`metadata`](@ref) function.

```@example 1
using NMRTools
# load an example 2D 1H,15N HMQC spectrum
spec = exampledata("2D_HN")

metadata(spec)
```

Entries can be accessed by passing a second argument to the `metadata` function:

```@example 1
metadata(spec, :title)
```

or directly from the spectrum using a dictionary-style lookup:

```@example 1
spec[:ns]
```

## Accessing axis metadata

Axis metadata can be accessed by providing an additional argument to the [`metadata`](@ref) function, specifying the axis numerically:

```@example 1
metadata(spec, 1)
```

Axes can also be accessed by their type, e.g. `F1Dim` or `F2Dim`:

```@example 1
metadata(spec, F2Dim)
```

Individual entries can be accessed by passing an additional argument to the `metadata` function:

```@example 1
metadata(spec, F1Dim, :label)
```

or directly from the spectrum object using a dictionary-style lookup alongside an axis reference:

```@example 1
spec[2, :offsetppm]
```


## Accessing acquisition parameters

Spectrometer acquisition parameters are automatically parsed from the `acqus` file when data are loaded. This is stored as a dictionary in the `:acqus` entry of the spectrum metadata, but can also be accessed through the function `acqus(spec, parametername)`. Parameter names can be provided either as strings (case insensitive, e.g. `"TE"` for the temperature) or as lowercase symbols (e.g. `:te`).

```@example 1
acqus(spec, :te)
```

Arrayed parameters such as pulse lengths are returned as dictionaries:

```@example 1
acqus(spec, :p)
```

Particular entries can be accessed directly by supplying an additional index parameter:

```@example 1
acqus(spec, :cnst, 4)
```

!!! note
    Arrayed parameters such as pulse lengths, delays, etc. are returned as dictionaries rather than lists to avoid indexing confusion between Bruker arrays, which are zero-based, and Julia arrays, which are unit-based.


## Spectrometer channels

Frequency dimensions carry their own `:bf`, `:sf`, and `:nucleus` metadata, but these only describe nuclei that were *detected* (i.e. that have an axis). To reach a nucleus that was manipulated but never detected — for example the carbon channel in a ¹³C-edited experiment read out on ¹H — use the spectrometer **channel** model.

Channels describe the physical RF channels (`:f1`…`:f8`) independently of the detected dimensions. `channels(spec)` returns a dictionary of all channels, and `channels(spec, ...)` returns a single channel by label, by `Nucleus`, or by a string giving either:

```@example 1
channels(spec, :f1)        # by spectrometer channel label
```

```@example 1
channels(spec, "13C")[:bf] # carbon base frequency, even with no carbon axis
```

Each channel dictionary has keys `:label`, `:nucleus`, `:bf` (base frequency, Hz), `:sfo` (carrier frequency, Hz), `:offsethz`, and `:offsetppm`.

## Auxiliary files

If present, files such as `vclist` and `vdlist` are imported and can be accessed through the `acqus` function:

```@example 1
relaxation_experiment = exampledata("pseudo3D_HN_R2")
acqus(relaxation_experiment, :vclist)
```


## Window functions

The window functions used for data processing are identified when experiments are loaded. These are represented by subtypes of [`WindowFunction`](@ref), which contain appropriate parameters specifying the particular function applied together with the acquisition time $t_\mathrm{max}$ (calculated at the point of application, i.e. after linear prediction but before zero filling).

Window functions are stored as axis metadata and can be accessed through the `:window` parameter:
```@example 1
spec[2, :window] # get the window function for the second dimension
```

Available window functions are:
- `NullWindow(tmax)`: no window function applied
- `UnknownWindow(tmax)`: unrecognised window function
- `ExponentialWindow(lb, tmax)`: line broadening specified in Hz
- `SineWindow(offset, endpoint, power, tmax)`: a supertype encompassing some special cases
    - `CosWindow(tmax)`: apodization by a pure half-cosine window
    - `Cos2Window(tmax)`: apodization by a pure half-cosine-squared window
    - `GeneralSineWindow(offset, endpoint, power, tmax)`: other general cases
- `GaussWindow(expHz, gaussHz, center, tmax)`: a supertype encompassing some special cases
    - `LorentzToGaussWindow(expHz, gaussHz, tmax)`
    - `GeneralGaussWindow(expHz, gaussHz, center, tmax)`: other general cases
