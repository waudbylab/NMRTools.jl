# Metadata

[`NMRData`](@ref) objects can store various metadata associated with the spectrum and each of the dimensions.

Metadata are stored as dictionaries using symbols as keys (e.g. `:ns`). They can be accessed using the `metadata` function, or directly from an `NMRData` object using a dictionary-style lookup. Metadata associated with axes are accessed by providing an additional reference, either as a dimension number or type (e.g. `F1Dim`, `F2Dim` etc.).

```julia
metadata(nmrdata, key) # spectrum metadata
nmrdata[key]

metadata(nmrdata, dimension, key) # axis metadata
nmrdata[dimension, key]
```

```@docs; canonical=false
metadata
```

## Labels

```@docs; canonical=false
label
label!
units
```

## Acquisition data

When spectra are loaded, the contents of the `acqus` file  are parsed (as are the `acqu2s` files etc. too, if present). These can be accessed with `:acqus` and `:acqu2s` keys, etc. For convenience though, the additional function [`acqus`](@ref) is provided to access acquisition data directly.

```@docs; canonical=false
acqus
```

## Auxiliary files

`acqus` can also be used to access the contents of auxiliary files (and if not present, `nothing` will be returned). Note that NMRTools will perform automatic unit conversion as follows:
- `:vclist`: variable loopcounter
- `:vdlist`: variable delays, in seconds
- `:valist`: variable amplitude, in dB (converted from Watts if necessary)
- `:vplist`: variable pulse lengths, in seconds
- `:fq1list` up to `:fq8list`: frequency lists – see the [Frequency lists](../metadata/fqlists.md) page for more information.

## Frequency list API

Frequency lists can be specified on the spectrometer in a number of ways - in Hz, in ppm, and relative to the spectrometer frequency or the base frequency (0 ppm). Frequency lists are therefore stored in NMRTools as `FQList` structures which encode this additional information.

```@docs; canonical=false
FQList
ppm
hz
```

Raw numerical data can be accessed using the `data()` function, but it is recommended to use `ppm` and `hz` functions to access frequency data safely. See the [Frequency lists](../metadata/fqlists.md) page for detailed examples.


## Standard metadata: spectra

| Key                  | Description                                                         |
|:-------------------- |:------------------------------------------------------------------- |
| `:acqusfilename`     | path to acqus file                                                  |
| `:acqus`             | contents of acqus file, as a dictionary                             |
| `:acqu2s`, `:acqu3s` | contents of acqu2s/acqu3s files (if present)                        |
| `:date`              | experiment timestamp (completion time)                              |
| `:experimentfolder`  | path to experiment                                                  |
| `:field`             | field strength, in T                                                |
| `:filename`          | spectrum filename or template                                       |
| `:format`            | input file format (`:nmrpipe`, `:pdata`, …)                         |
| `:label`             | short label (first line of title pdata file), used for captions     |
| `:ndim`              | number of dimensions                                                |
| `:noise`             | RMS noise level (see [`estimatenoise`](@ref))                       |
| `:nuclei`            | set of `Nucleus` values                                             |
| `:ns`                | number of scans                                                     |
| `:pulseprogram`      | pulse programme (PULPROG)                                           |
| `:rg`                | receiver gain                                                       |
| `:solvent`           | solvent string, e.g. `"H2O+D2O"`                                    |
| `:temperature`       | temperature, in K                                                   |
| `:title`             | spectrum title (contents of title pdata file)                       |
| `:topspin`           | TopSpin version used for acquisition                                |

## Standard metadata: frequency and time dimensions

Keys marked † are set on time dimensions by [`loadfid`](@ref NMRTools.NMRIO.loadfid) and carried
over to the frequency dimension by [`ft`](@ref).

| Key                  | Description                                                                |
|:-------------------- |:-------------------------------------------------------------------------- |
| `:aq`                | acquisition time, in seconds                                               |
| `:bf`                | base frequency, in Hz                                                     |
| `:firstpoint` †      | scaling factor applied to the t = 0 point                                  |
| `:grpdly` †          | digital filter group delay, in dwell times (direct dimension; 0 after `ft`) |
| `:label`             | short label                                                                |
| `:mcindex`           | imaginary unit of the dimension (`1` for `im1`, …), `nothing` if none      |
| `:npoints`           | final number of **real** data points in dimension (after extraction)       |
| `:nucleus`           | `Nucleus` of the dimension                                                 |
| `:offsethz`          | carrier offset from bf, in Hz                                              |
| `:offsetppm`         | carrier offset from bf, in ppm                                             |
| `:phase`             | total phase correction applied, `(p0=…, p1=…)` in degrees, p0 at the first point |
| `:pseudodim`         | flag indicating non-frequency domain data (`false` for frequency domain)   |
| `:quadrature` †      | quadrature scheme of an indirect dimension, e.g. `:states_tppi`            |
| `:referenceoffset`   | applied referencing, in ppm                                                |
| `:region`            | extracted region, expressed as a range in points, otherwise `missing`      |
| `:sf`                | carrier frequency, in Hz                                                  |
| `:swhz`              | spectrum width, in Hz                                                      |
| `:swppm`             | spectrum width, in ppm                                                     |
| `:td`                | number of **complex** points acquired, including LP                        |
| `:tdzf`              | number of **complex** points when FT executed, including LP and ZF         |
| `:window`            | [`WindowFunction`](@ref) encoding applied apodization                      |
