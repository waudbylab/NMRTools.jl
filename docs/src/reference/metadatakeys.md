# Metadata keys

The standard metadata keys set when data are loaded and processed. See [Metadata](../concepts/metadata.md)
for how metadata are organised.

## Accessing metadata

| Lookup | Returns |
|:-------|:--------|
| `spec[:ns]` | spectrum metadata (a key from [Spectra](#Spectra)) |
| `spec[F1Dim, :swppm]`, `spec[1, :swppm]` | axis metadata (a key from [Frequency and time dimensions](#Frequency-and-time-dimensions)), by dimension type or number |
| `acqus(spec, :te)` | an acquisition parameter from the `acqus` file |
| `acqus(spec, :p, 1)` | an arrayed acquisition parameter (here `P1`), indexed as in TopSpin |
| `acqus(spec, :vdlist)` | an auxiliary file (see [Auxiliary files](#Auxiliary-files)) |
| `metadata(spec)`, `metadata(spec, F1Dim)` | the whole dictionary for the spectrum, or for one axis |

```@example keys
using NMRTools
spec = exampledata("2D_HN")
spec[:ns], spec[F1Dim, :swppm], acqus(spec, :p, 1)
```

Times are in seconds and frequencies in Hz (see [Time and frequency units](../concepts/units.md)).

## Spectra

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

## Frequency and time dimensions

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


## Auxiliary files

Auxiliary files in the experiment folder are read when data are loaded and returned by
[`acqus`](@ref), with their values converted to NMRTools units:

| Key | Contents |
|:--- |:-------- |
| `:vclist` | variable loop counters (dimensionless) |
| `:vdlist` | variable delays, in seconds |
| `:valist` | variable amplitudes, in dB (converted from W if necessary) |
| `:vplist` | variable pulse lengths, in seconds |
| `:fq1list` … `:fq8list` | frequency lists, as [`FQList`](@ref) (see [Frequency lists](fqlists.md)) |

`acqus` returns `nothing` for a list that is not present.
