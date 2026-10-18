# Working with NMR data

NMR measurements are arrays of data, with additional numerical data associated with each dimension, or axis. Within NMRTools, these data are stored as [`NMRData`](@ref) structures, which encapsulate the data, axis information, and metadata on acquisition or processing.


## Loading NMR data

NMR data are loaded using the [`loadnmr`](@ref) function. This handles processed Bruker experiments and NMRPipe-format data.

```julia
# load bruker experiment number 1 from a directory '2D_HN'
# by default, NMRTools will load bruker processed data from proc 1
spec2d = loadnmr("exampledata/2D_HN/1")

# load a different processed spectrum
spec1d = loadnmr("exampledata/2D_HN/1/pdata/101")

# load data from NMRPipe format, using a template
spec3d = loadnmr("exampledata/pseudo3D_HN_R2/1/ft/test%03d.ft2")
```

!!! tip
    `loadnmr` will attempt to locate and parse acquisition metadata, such as acqus files. If the spectrum file is located elsewhere (for example, if you are loading a file that was processed with NMRPipe), then you can specify the path to the experiment folder using the `experimentfolder` keyword argument.

When spectra are loaded, a simple algorithm runs to estimate the noise level, which is often used for subsequent plotting commands.

## Manipulating spectrum data

[`NMRData`](@ref) structures encapsulate a standard Julia array. This can be accessed using the [`data`](@ref) command. Through multiple dispatch, most operations work directly on NMRData variables as if they are regular arrays, with the added benefit that axis information and metadata are preserved. Data can be sliced and accessed like a regular array using the usual square brackets:
```julia
spec1d[100:105]
spec2d[3:4, 10:14]
```

Value-based selectors can also be used to locate data using chemical shifts. Three selectors are defined:
- `At(x)`: select data precisely at the specified value
- `Near(x)`: select data at the nearest matching position
- `x .. y`: select the range of data between `x` and `y` (closed interval)

For example:
```julia
spec1d[8.2 .. 8.3] # select between 8.2 and 8.3 ppm
spec2d[Near(8.25), 123 .. 124] # select near 8.25 ppm in the first dimension
                               # and between 123 and 124 ppm in the second dimension
```

When data are sliced, new NMRData structures are created and their axes are updated to match the new data size.

!!! warning
    When `NMRData` structures are sliced, copied, or otherwise modified, they inherit the same dictionary of metadata as the original variable. This means that any changes to metadata will affect both variables. To resolve this, make a `deepcopy` of the variable. Note also that any acquisition metadata might not reflect the correct shape of the data any more.


## Accessing multicomplex (real/imaginary) data

NMR data can be loaded with all quadrature components using the `allcomponents=true` keyword:

```julia
spec1d = loadnmr("expt/110", allcomponents=true)
spec2d = loadnmr("expt/111", allcomponents=true)
```

Complex multidimensional data is represented using [MulticomplexNumbers.jl](https://waudbylab.org/MulticomplexNumbers.jl). This algebra provides a natural representation of the multiple real and imaginary components associated with NMR data ([Delsuc 1988](https://www.sciencedirect.com/science/article/pii/0022236488900364)), and facilitates phase correction and bookkeeping of individual real and imaginary components.

A 2D point has four components, `rr + ir·im1 + ri·im2 + ii·im1·im2`, where `im1` belongs to the direct dimension and `im2` to the indirect one. They can be extracted with `component`, numbered in that order; `realest(spec2d)` gives the purely real part `rr`:

```julia
fig = Figure(size=(600, 600))
nmrplot(fig[1,1:2], [real(spec1d), imag(spec1d)])  # 1D plot
spec2d = spec2d[-0.5 .. 1.2, 81 .. 90] / 2
nmrplot(fig[2,1], component.(spec2d, 1), title="component 1 (rr)")  # 2D plots
nmrplot(fig[2,2], component.(spec2d, 2), title="component 2 (ir)")  # imaginary in the direct dimension
nmrplot(fig[3,1], component.(spec2d, 3), title="component 3 (ri)")  # imaginary in the indirect dimension
nmrplot(fig[3,2], component.(spec2d, 4), title="component 4 (ii)")
```

Imaginary components are the complex conjugates of TopSpin's (`2ir`, `2ri`, …), following the NMRTools sign convention.

![Multicomplex data](../assets/multicomplex.png)

## Accessing axis data

Information on data dimensions is stored in `NMRDimension` structures. These can be accessed with the `dims` function:
```@example 1
# get the first dimension of this two-dimensional experiment
using NMRTools # hide
spec2d = exampledata("2D_HN"); # hide
dims(spec2d, 1)
```

`NMRDimension`s can be treated like vectors (one-dimensional arrays) for most purposes, including indexing and slicing. Value-based selectors can also be used, as for spectrum data. Like spectrum data, the underlying numerical data can be accessed if needed using the [`data`](@ref) function.

A hierarchy of types is defined for NMR dimensions, reflecting the variety of different experiments:
- `NMRDimension`
    - `FrequencyDimension`: with specific types `F1Dim` to `F4Dim`
    - `NonFrequencyDimension`
        - `TimeDimension`:
            - `TrelaxDim`: for relaxation times
            - `TkinDim`: for kinetic evolution times
            - `T1Dim` to `T4Dim`: for general frequency evolution periods
        - `GradientDimension`: for e.g. diffusion measurements, with specific types `G1Dim` to `G4Dim`
        - `UnknownDimension`: with specific types `X1Dim` to `X4Dim`

## Chemical shift referencing

The [`reference`](@ref) function corrects the referencing of a spectrum by specifying a known chemical shift: you provide an `old => new` pair indicating where a reference peak currently appears and where it should appear.

```@example offset
using NMRTools, Plots #hide

# Load example data
spec2d_original = exampledata("2D_HN")

# Reference the 1H dimension: peak at 4.72 ppm should be at 4.70 ppm
spec2d = reference(spec2d_original, F2Dim, 4.72 => 4.70)

label!(spec2d_original, "unreferenced")
label!(spec2d, "referenced")
plot([spec2d_original, spec2d])
savefig("plot-offset.svg"); nothing # hide
```

![](plot-offset.svg)

The dimension can be specified as a dimension index (`1`, `2`, ...), a dimension type (`F1Dim`, `F2Dim`, ...), or a nucleus (`H1`, `N15`, ...).

By default, `reference` also applies indirect referencing to the other frequency dimensions using IUPAC Xi ratios, so correcting the ¹H reference of an HSQC corrects its ¹⁵N or ¹³C axis too. For aqueous samples, `reference(spec)` references to the water signal using the sample temperature. See [Chemical Shift Referencing](@ref) for these options and the Xi ratio tables.


## Accessing metadata

[`NMRData`](@ref) objects contain metadata on processing and acquisition parameters that are populated automatically upon loading a spectrum. Entries are divided into *spectrum* metadata - associated with the experiment in general - and *axis* metadata, that are associated with a particular dimension of the data.

Metadata entries are labelled by symbols such as `:ns` or `:pulseprogram`. Entries can be accessed using the `metadata` function, or directly as a dictionary-style lookup:
```@repl 1
spec2d = exampledata("2D_HN"); # hide
metadata(spec2d, :ns)
spec2d[:title]
```

Acquisition parameters from Bruker acqus files are also parsed when loading data, and can be accessed using the [`acqus`](@ref) function. Parameters are specified either as lower-case symbols or strings (not case-sensitive). An index can be specified for arrayed parameters such as pulse lengths or delays.
```@repl 1
acqus(spec2d, "TE")
acqus(spec2d, :p, 1)
```

Axis metadata can be accessed by providing an additional label, which can either be numerical or the appropriate `NMRDimension` type, such as `F1Dim` etc:
```@repl 1
metadata(spec2d, F2Dim, :label)
spec2d[2, :bf]
spec2d[F2Dim, :window]
```

See the [Metadata overview](@ref) for more on metadata, including spectrometer channels, power levels, frequency lists and sample information.


## Saving data

Frequency-domain spectra can be written as Bruker processed data with [`savenmr`](@ref), and time-domain data as Bruker `fid` or `ser` files with [`savefid`](@ref), for example to view a spectrum processed in NMRTools in TopSpin:

```julia
savenmr("myexperiment/1", spec; procno=101)   # writes myexperiment/1/pdata/101
```

Peak lists are saved with [`savepeaks`](@ref) (see [Peaks and peak lists](peaks.md)).
