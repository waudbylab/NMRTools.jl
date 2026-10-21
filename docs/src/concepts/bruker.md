# Bruker data and TopSpin conventions

NMRTools reads and writes Bruker raw and processed data, and is designed to give the same
spectra as TopSpin. This page describes what [`loadfid`](@ref NMRTools.NMRIO.loadfid) does to raw
data, how the digital filter is handled, and the conventions that differ between the two programs.

## What `loadfid` does

On loading a `fid` or `ser` file:

- indirect dimensions acquired with States, States-TPPI or echo/anti-echo are recombined into an
  extra multicomplex order, so a 2D experiment has `Multicomplex{Float64,2,4}` elements; the scheme
  used is recorded in each dimension's `:quadrature` metadata;
- non-phase-sensitive dimensions (relaxation delays, CEST offsets, …) become `X2Dim`/`X3Dim`;
- the digital filter's build-up points stay in the FID, on a time axis that starts before t = 0
  (see below);
- the t = 0 point of each Fourier dimension is scaled by 0.5 (change this with the `firstpoint`
  keyword), except in the direct dimension of baseopt data;
- the data are scaled by `2^NC`, as TopSpin scales them, so that intensities after [`ft`](@ref),
  which applies no normalisation, match TopSpin's.

[`savefid`](@ref NMRTools.NMRIO.savefid) undoes each of these steps, so `loadfid` on a written
experiment returns the same data.

## Digital filtering and group delay

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

### Baseopt

In baseopt mode the spectrometer aligns acquisition with time zero (the pulse programme's `acqt0`)
and corrects the initial points of the FID itself, so no first-order phase correction should be
needed and the t = 0 point already has its half weight.

### FIDs modified by `cryoproc1d`

A `DSPFVS` outside 10–23 is not a firmware version Bruker defines. TopSpin's `cryoproc1d`, used to
process X-nucleus spectra from cryoprobes, sets `DSPFVS = 0` after rewriting the FID with its start
removed: the filter build-up and some points after time zero (the original is saved in another
experiment number, recorded in `audita.txt`). `loadfid` compares `TD` in `acqu`, which keeps the
original value, with `TD` in `acqus` to find how many points after time zero are missing, and starts
the time axis at that time, so that `ft` restores the correct time origin and the lines can be
phased together. A warning reports the number of missing points.

Their absence leaves a ripple around strong lines. TopSpin rebuilds them by backward linear
prediction with thousands of coefficients; `linearpredict(fid, T1Dim; direction=:backward,
npoints=n)`, or `autoprocess(fid; rebuild=true)`, rebuilds them with a low-order prediction, which
removes the ripple but can distort the baseline when the FID contains fast-decaying signals such as
probe background.

## Phase conventions

`phase(spec, dim, p0, p1)` applies the phase `θ = p0 + p1·n/N` to point `n` (counting from zero at
the left-hand, highest-ppm edge), in degrees. The sign and pivot match TopSpin's `PHC0` and `PHC1`,
so values can be copied between the two. Each call adds to the dimension's `:phase` metadata, which
holds the total correction applied so far as `(p0=…, p1=…)` with `p0` referred to the first point.

!!! note "Echo/anti-echo data"
    For 1D data, NMRTools and TopSpin need the same `PHC0`/`PHC1`. For 2D echo/anti-echo data the
    phases currently differ from TopSpin's, by about 90° in the direct dimension and 180° in the
    indirect dimension. This is a known issue in the echo/anti-echo recombination in `loadfid`.

## Sign convention

NMRTools spectra are the complex conjugate of TopSpin's in every imaginary unit: the real parts
agree and the imaginary parts have the opposite sign. [`phase`](@ref) uses TopSpin's sign for the
phase angle, which is why `PHC0` and `PHC1` carry over unchanged.

When imaginary data are read from or written to Bruker processed files, they are converted between
the two conventions: `loadnmr(...; allcomponents=true)` conjugates the imaginary components it
reads, and [`savenmr`](@ref NMRTools.NMRIO.savenmr) conjugates those it writes. Each program's
files therefore follow its own convention. For 2D data, `2ir` holds the coefficient of `im1` (the
direct dimension) and `2ri` that of `im2`.

## Intensities

`loadfid` scales raw data by `2^NC` and `ft` applies no normalisation, so spectra processed in
NMRTools have the same intensities as TopSpin's (to about 1 part in 10⁸ for the `1r` of a 1D).
`savenmr` writes 32-bit integer data and sets `NC_proc` so that TopSpin reads the same intensities.
