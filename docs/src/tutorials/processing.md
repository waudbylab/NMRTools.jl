# Processing raw data

!!! warning "Under development"
    Loading and processing raw time-domain data is new and its interface may change. The examples on
    this page use raw Bruker data, which is not shipped with the package, so they are not run when
    the documentation is built.

NMRTools can process raw Bruker data (`fid`/`ser` files) into spectra. Processing steps are ordinary
functions that take and return [`NMRData`](@ref), so they can be chained with `|>`. Throughout
processing the element type is `Multicomplex`: each dimension's quadrature is carried by its own
imaginary unit (`im1` for the direct dimension, `im2` for the first indirect dimension, …), and the
`:mcindex` metadata of each dimension records which unit belongs to it.

## Quick start: `autoprocess`

[`autoprocess`](@ref) turns raw data into a phased spectrum in one step: water suppression for ¹H
spectra in H₂O, a cos² window, zero-filling to twice the size, Fourier transformation and automatic
phasing of every dimension.

```julia
using NMRTools
spec = loadfid("myexperiment/1") |> autoprocess
spec = loadfid("myexperiment/1") |> autoprocess(window=ExponentialWindow(1.0), solvent=false)
```

Each step can be changed with a keyword, for all dimensions or per dimension with a dictionary:

```julia
spec = autoprocess(loadfid("hsqc/1");
                   window=Dict(T1Dim => Cos2Window(), T2Dim => Cos2Window()),
                   zerofill=Dict(T2Dim => 4),
                   phase=Dict(F2Dim => (0.0, 0.0)))     # fixed phases for the indirect dimension
```

The rest of this page describes the individual steps.

## Loading raw data

[`loadfid`](@ref) reads a `fid` (1D) or `ser` (nD) file and returns time-domain data:

```julia
using NMRTools
fid = loadfid("myexperiment/1")       # or the path of the fid/ser file
```

On loading:

- indirect dimensions acquired with States, States-TPPI or echo/anti-echo are combined into an extra
  multicomplex order, so a 2D experiment has `Multicomplex{Float64,2,4}` elements;
- non-phase-sensitive dimensions (relaxation delays, CEST offsets, …) become `X2Dim`/`X3Dim`;
- the digital filter's build-up points stay in the FID, on a time axis that starts before t = 0;
  [`ft`](@ref) removes the group delay (see below);
- the t = 0 point of each Fourier dimension is scaled by 0.5 (change this with `firstpoint`), except
  in the direct dimension of baseopt data;
- the data are scaled by `2^NC`, as TopSpin scales them, so that intensities after [`ft`](@ref),
  which applies no normalisation, match TopSpin's.

### Digital filtering and baseopt

Bruker spectrometers filter and decimate the oversampled signal digitally. The filter's response
builds up over the first `GRPDLY` points of the stored FID, so time zero falls `GRPDLY` points into
the data. The filter's response is symmetric, so those build-up points carry part of the signal at
negative time.

`loadfid` keeps the build-up points: the direct dimension's time axis starts at `-GRPDLY` dwell
times, and the delay is recorded in its `:grpdly` metadata. The processing functions respect this:

- [`apodize`](@ref) starts the window at t = 0 and leaves the build-up points unwindowed;
- [`zerofill`](@ref) appends zeros after the end of the FID;
- [`linearpredict`](@ref) predicts from the FID from t = 0;
- [`ft`](@ref) removes the delay with a linear phase ramp after the Fourier transform, which moves
  the build-up points to negative time (the end of the zero-filled array).

This is how TopSpin and NMRPipe process the data, and it gives the same spectrum as TopSpin.
Discarding the build-up points instead truncates the filter's response and produces a "smile": a
baseline that falls away at both edges of the spectrum.

How the delay is found depends on the acquisition parameters:

| Parameters | Treatment |
|---|---|
| `DIGMOD = 3` (baseopt, `DSPFIRM` rectangle) | delay `GRPDLY` (an integer); t = 0 point not scaled |
| `DIGMOD = 1` (digital), `DSPFVS` 20–23 | delay `GRPDLY` (may be fractional); t = 0 point scaled by 0.5 |
| `DIGMOD = 0` (analogue) | no delay |
| `DSPFVS` outside 10–23 (FID modified, for example by `cryoproc1d`) | time axis starts after t = 0 by the points removed beyond `GRPDLY`, with a warning |

In baseopt mode the spectrometer aligns acquisition with time zero (the pulse programme's `acqt0`)
and corrects the initial points of the FID itself, so no first-order phase correction should be
needed and the t = 0 point already has its half weight.

A `DSPFVS` outside 10–23 is not a firmware version Bruker defines. TopSpin's `cryoproc1d`, used to
process X-nucleus spectra from cryoprobes, sets `DSPFVS = 0` after rewriting the FID with its start
removed: the filter build-up and some points after time zero (the original is saved in another
experiment number, recorded in `audita.txt`). `loadfid` compares `TD` in `acqu`, which keeps the
original value, with `TD` in `acqus` to find how many points after time zero are missing, and starts
the time axis at that time, so that `ft` restores the correct time origin and the lines can be
phased together. A warning reports the number of missing points. Their absence leaves a ripple
around strong lines. TopSpin rebuilds them by backward linear prediction with thousands of
coefficients; `linearpredict(fid, T1Dim; direction=:backward, npoints=n)`, or
`autoprocess(fid; rebuild=true)`, rebuilds them with a low-order prediction, which removes the ripple
but can distort the baseline when the FID contains fast-decaying signals such as probe background.

## The processing pipeline

```julia
spec = loadfid("myexperiment/1") |>
       apodize(ExponentialWindow(1.0)) |>
       zerofill |>
       ft |>
       phase(F1Dim, 128.5) |>
       real
```

| Function | Acts on | Effect |
|---|---|---|
| [`apodize`](@ref) | time domain | multiply by a window function |
| [`zerofill`](@ref) | time domain | append zeros (default: double the size) |
| [`linearpredict`](@ref) | time domain | extend a dimension, or repair its first points |
| [`solventfilter`](@ref) | time domain | remove a solvent signal |
| [`ft`](@ref) | either | Fourier transform (`inverse=true` for the inverse) |
| [`phase`](@ref) | either | zero- and first-order phase correction |
| [`autophase`](@ref) | frequency domain | estimate and apply a phase correction |
| `real` | either | take the real part, giving `Float64` data |

With no dimension argument, a function acts on every Fourier dimension. Pass a dimension to act on
one only, for example to use different windows in each dimension of a 2D experiment:

```julia
spec = loadfid("hsqc/1") |>
       apodize(T1Dim, Cos2Window()) |> apodize(T2Dim, Cos2Window()) |>
       zerofill |> ft |>
       phase(F1Dim, -87.5) |> phase(F2Dim, 0.0, 180.0) |>
       real
```

Each function also accepts a `Dict` of spectra and processes every entry.

### Phase conventions

`phase(spec, dim, p0, p1)` applies the phase `θ = p0 + p1·n/N` to point `n` (counting from zero at
the left-hand, highest-ppm edge), in degrees. The sign and pivot match TopSpin's `PHC0` and `PHC1`,
so values can be copied between the two. Pass `pivot` (in ppm) to move the point at which only `p0`
acts.

Each call adds to the dimension's `:phase` metadata, which holds the total correction applied so far
as `(p0=…, p1=…)` with `p0` referred to the first point:

```julia
metadata(spec, F1Dim, :phase)    # (p0 = 128.5, p1 = 0.0)
```

!!! note "Echo/anti-echo data"
    For 1D data, NMRTools and TopSpin need the same `PHC0`/`PHC1`. For 2D echo/anti-echo data the
    phases currently differ from TopSpin's, by about 90° in the direct dimension and 180° in the
    indirect dimension. This is a known issue in the echo/anti-echo recombination in `loadfid`.

## Automatic phasing

[`autophase`](@ref) estimates the phase correction of a frequency-domain spectrum and applies it:

```julia
spec = loadfid("myexperiment/1") |> apodize(ExponentialWindow(1.0)) |> zerofill |> ft |>
       autophase |> real
# [ Info: autophase: p0 = 141.3°, p1 = -4.2° (residual 11.0° over 15 line groups)
```

To see the estimate without applying it, use [`estimatephase`](@ref). It returns the phases with a
residual (the weighted rms phase error of the lines used, in degrees), the number of line groups
used, and the groups left out of the fit because their phase disagreed. A large residual, few
groups, or a strong rejected group indicates an unreliable estimate:

```julia
est = estimatephase(spec)
est.p0, est.p1, est.residual, est.nlines
est.rejected      # [(shift = -130.36, residual = 70.0, relsnr = 0.029)]
est.edges         # [(edge = :left, snr = 95.2, residual = 0.4, used = true), …]
phase(spec, F1Dim, est.p0, est.p1)
```

`autophase` warns when a rejected line is at least 2% as strong as the strongest line used. Such a
line may have a genuinely different phase, or the first-order phase may be far from its expected
value. The second case is common for heteronuclei with wide spectral windows: there, set the
expected first-order phase with `p1` (increasing `maxdelay` if needed) or fix it with `order=0`.

Useful options:

- `order=0` fits the zero-order phase only, with the first-order phase fixed at `p1` (default 0).
  Once the digital-filter group delay has been removed, the first-order phase of the direct
  dimension should be close to zero, so this is often the more robust choice.
- `maxdelay` sets how far the first-order phase may range. A first-order phase is a shift of the
  time origin, so the range is given as a timing uncertainty (default 20 µs) and converted with the
  spectral width: ±100° for a 14 kHz ¹H spectrum, ±950° for a 130 kHz ¹⁹F spectrum.
- `exclude=[4.6 .. 4.9]` ignores lines in the given chemical-shift ranges, for example water or buffer.
- `λ` sets how strongly the first-order phase is pulled towards `p1`. The pull is scaled by the
  scatter of the line phases and the number of lines, so it has little effect when many lines agree.
  Increase it when the lines cover only part of the spectrum.

### How it works

The spectrum is modelled as a sum of Lorentzian lines with [`estimatepoles`](@ref) (see below), so
each line has a frequency and a phase. Lines within a coupling constant of each other are grouped
and the phase of their summed amplitude is used: echoes in the pulse sequence (excitation
sculpting, for example) give the lines of a multiplet antiphase dispersive components that cancel
only in the sum. Lines much broader than the typical linewidth are ignored when they overlap a
narrower line, since they then describe baseline distortions or the tails of a non-Lorentzian line
such as a radiation-damped water signal. Isolated broad lines are kept.

The two edges of the spectrum give two further observations. With the correct phase the real part
of the spectrum is flat far from the lines, while the imaginary part carries the dispersion tails of
the strong lines. A small phase error leaks those tails into the real part, tilting the baseline at
the edges (zero order) or curling it (first order). The phase at which the real part of each edge
region vanishes is therefore a precise phase measurement wherever the edges carry signal, as they do
in protein spectra and spectra with residual water.

The correction is a weighted circular regression of the observed phases against position. Each
observation is weighted by its inverse variance: a statistical part from its signal-to-noise ratio,
plus a systematic part, which for line groups is their scatter about the fit and for edges is 2°.
The line phases are taken modulo 180°, so negative peaks are treated as correctly phased rather than
as phase errors. For each first-order phase the best zero-order phase has a closed form, so the fit
is a one-dimensional search over the first-order phase. Observations that disagree with the fit by
more than `maxresidual` degrees are discarded and the fit repeated. The remaining 180° ambiguity is
resolved by a majority vote of the line groups, so a single dominant line (residual water, say)
cannot decide it.

Compared with hand phasing in TopSpin on ¹H spectra of proteins and small molecules, the estimates
agreed within a few degrees across the spectrum where the edges carry signal. Without edge signal,
as in a broadband ¹⁹F spectrum with weak lines, the estimate rests on the lines alone and can be off
by 10–30° at the edges. With only one or two lines, the first-order phase cannot be determined from
the lines and stays close to `p1`. Check the result on new types of data, and use `order=0`, `p1`
and `exclude` where the fit is poor.

For multidimensional data, `autophase(spec, dim)` estimates the correction from a single trace along
`dim` and applies it to the whole dataset. The trace takes the first point of every other time
dimension and sums over the remaining dimensions.

## Pole models

[`estimatepoles`](@ref) models a spectrum as a sum of damped complex exponentials, one
[`Pole`](@ref) per line, using band-wise harmonic inversion (HSVD). Each pole has a chemical shift,
a linewidth (Hz), a complex amplitude (whose magnitude is the line integral and whose argument is
the line's phase) and a signal-to-noise ratio:

```julia
poles = estimatepoles(spec)
# Pole(3.6848 ppm, 2.6 Hz, phase 128.2°, SNR 2.31e4)
```

## Linear prediction and solvent suppression

[`linearpredict`](@ref) extends a time-domain dimension by linear prediction, typically a truncated
indirect dimension. Process the direct dimension first, so that each indirect vector belongs to one
direct-dimension frequency:

```julia
spec = loadfid("hsqc/1") |>
       apodize(T1Dim) |> zerofill(T1Dim) |> ft(T1Dim) |>
       linearpredict(T2Dim) |>
       apodize(T2Dim) |> zerofill(T2Dim) |> ft(T2Dim) |>
       real
```

`linearpredict(spec, dim; replace=n)` replaces the first `n` points by backward prediction instead.

[`solventfilter`](@ref) removes a solvent signal from the FID before apodization. By default it acts
on the carrier frequency; pass `shift` (ppm) to target another position. `method=:lowpass` subtracts
a low-pass filtered copy of the FID (the convolution difference method of nmrPipe's `SOL`), and
`method=:hsvd` fits and subtracts only the lines within a band around the target, which is slower but
more selective.

## Saving in Bruker format

[`savenmr`](@ref) writes a frequency-domain spectrum as Bruker processed data, and [`savefid`](@ref)
writes time-domain data as a raw Bruker `fid`/`ser`:

```julia
fid = loadfid("myexperiment/1")
spec = fid |> apodize(ExponentialWindow(1.0)) |> zerofill |> ft |> autophase

savenmr("myexperiment/1", spec; procno=2)          # writes myexperiment/1/pdata/2
savefid("processed/1", linearpredict(fid))          # new experiment folder
```

`savenmr` writes the real part (`1r`, `2rr`, `3rrr`) and, for multicomplex 1D and 2D data, the
imaginary parts, as 32-bit integers scaled by `NC_proc`. The processing parameter files (`procs`,
`proc2s`, …) record the axes, the window function and the phase correction applied.

`savefid` undoes the steps taken by `loadfid`: it restores the first points and the digital-filter
group delay, re-encodes each indirect dimension with its original quadrature scheme and puts the
dimensions back in acquisition order. `loadfid` on the written experiment returns the same data.
Data are written as `Float64` (`DTYPA=2`, as recorded by Avance NEO spectrometers) by default; pass
`dtype=Int32` for older TopSpin versions.

Both functions copy the acquisition parameter files (`acqus`, `acqu2s`, pulse programme, lists, …)
from the experiment the data were loaded from, updating only the parameters that describe the
written data, such as `TD`. If there is no source experiment, minimal parameter files are generated
from the dimension metadata. Existing data are not overwritten unless `overwrite=true` is passed.

!!! note "Sign convention"
    NMRTools spectra are the complex conjugate of TopSpin's: the real parts agree and the imaginary
    parts have the opposite sign. `savenmr` and `loadnmr(...; allcomponents=true)` convert between
    the two, so imaginary data read from or written to Bruker files follow each program's own
    convention.
