# Time and frequency units

NMRTools stores **all times in seconds and all frequencies in Hz**. Bruker parameter files use a
mixture of units (pulse lengths in µs, some delays in ms, spectrometer frequencies in MHz), so
values are converted when acquisition parameters are read. The rule applies to:

- pulse lengths (`P`)
- delays (`D`), vdlists and vplists
- acquisition times and dwell times
- spectrometer base and carrier frequencies (`:bf`, `:sf`)
- spectral widths (`:swhz`) and carrier offsets (`:offsethz`)
- linewidths of peaks

Chemical shifts are the exception: they are in ppm, and frequency dimensions have ppm coordinates.
Time dimensions have coordinates in seconds.

```@example units
using NMRTools
spec = exampledata("2D_HN")
acqus(spec, :p, 1)          # pulse length p1, in s
```

```@example units
acqus(spec, :d, 1)          # delay d1, in s
```

```@example units
spec[F1Dim, :bf]            # ¹H base frequency, in Hz
```

## Converting between ppm and Hz

`ppm` and `hz` convert a frequency list to chemical shifts or to offsets from the carrier, using the
base frequency and carrier of a given dimension; see [Frequency lists](../reference/fqlists.md).

Power levels have their own type, [`Power`](@ref), which converts between dB and W; see
[Power levels](../reference/power.md).
