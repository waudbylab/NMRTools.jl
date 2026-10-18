# NMRTools.jl release notes

## v0.3.0

This release adds loading and processing of raw time-domain data, writing of Bruker data, a new
peak framework and import of assigned chemical shifts from the BMRB. It contains breaking changes,
listed first.

### Breaking changes

- **`loadnmr(...; allcomponents=true)` returns different numbers for the same files.** Imaginary
  components are now conjugated relative to TopSpin's files, and for 2D data the `2ir` and `2ri`
  files have swapped places: `2ir` is now the coefficient of `im1` (the direct dimension) and `2ri`
  that of `im2`. The element type is unchanged, so existing code will run without error but give
  different results. This follows the sign convention used throughout NMRTools: spectra are the
  complex conjugate of TopSpin's in every imaginary unit (real parts agree), so TopSpin's `PHC0` and
  `PHC1` can be passed straight to `phase`. Real data (the default, `allcomponents=false`) are
  unaffected.
- **Peaks.** `SimplePeak` and `AbstractPeak` are replaced by `Peak{N}`, holding a position in ppm,
  an intensity, a linewidth in Hz for each dimension and a label. `detectpeaks` returns a
  `Vector{Peak{N}}` and now works in any number of frequency dimensions.
- **`detectpeaks` keywords and defaults.** `snr_threshold` is renamed `snr` (the old name still
  works, with a deprecation warning, and will be removed in a future release). Negative peaks are
  now detected as well as positive ones by default; pass `sign=:positive` for the previous
  behaviour. Intensities are heights above the median baseline.
- **Dependencies.** Peaks.jl, ArtifactUtils.jl and VisualRegressionTests.jl are no longer
  dependencies of NMRTools.
- **Julia 1.10 or later** is now required (previously 1.9).
- **DimensionalData.jl 0.30** is now required (previously 0.29). DimensionalData 0.30 contains
  breaking changes of its own (see its
  [release notes](https://github.com/rafaqz/DimensionalData.jl/releases/tag/v0.30.0)), which
  affect code that uses DimensionalData functions directly on `NMRData`: for example
  `DimensionalData.LookupArrays` is now `DimensionalData.Lookups`, and the `rtol` keyword of `At`
  has been removed.

### Raw time-domain data

Loading and processing raw data is new, and its interface is the most likely part of the package to
change in a future release.

- `loadfid` reads Bruker `fid` and `ser` files (1D to 3D) as `Multicomplex` data, with one
  imaginary unit per Fourier dimension recorded in its `:mcindex` metadata. It combines States,
  States-TPPI and echo/anti-echo quadrature, scales the t = 0 point, and keeps the digital filter's
  build-up points on a time axis starting before t = 0. Intensities match TopSpin's. `loadnmr` on a
  `fid` or `ser` path now calls `loadfid`.
- Processing functions, each with a curried form for use with `|>`: `apodize`, `zerofill`, `ft`
  (including the inverse transform), `phase`, `extract`, `frequencyshift`, `linearpredict` (forward
  and backward) and `solventfilter`. `ft` removes the digital filter's group delay with a linear
  phase ramp, as TopSpin and NMRPipe do. `extract` selects a region, as `extract(spec, F1Dim,
  6.0 .. 10.0)` or in a pipeline as `extract(F1Dim, 6.0 .. 10.0)`.
- `real(spec)` gives the purely real part of a multicomplex spectrum (real in every dimension), and
  `imag(spec, dim)` the part imaginary in dimension `dim` alone (`imag(spec)` for 1D data);
  `abs` and `complex` are also defined.
- Automatic phasing: `estimatepoles`, `estimatephase`, `autophase` and `autophase0`. Without a
  dimension argument, `autophase` and `autophase0` phase every frequency dimension in turn (a
  dimension with too little signal is left unphased, with a warning). The total correction applied
  is recorded in each frequency dimension's `:phase` metadata, and the estimate is logged at debug
  level.
- `autoprocess` chains solvent suppression, apodization, zero-filling, Fourier transformation and
  automatic phasing in one step, with per-dimension settings.
- `Cos2Window` replaces `Cos²Window`, which remains as an alias. `CosWindow`, `Cos2Window` and
  `GeneralSineWindow` can be created without a length, which `apodize` then sets from the
  acquisition time.

### Writing data

- `savenmr(path, spec; procno=1)` writes frequency-domain data as Bruker processed data
  (`pdata/<procno>/1r`, `2rr`, …).
- `savefid(path, fid)` writes time-domain data as a Bruker `fid` or `ser` file, undoing everything
  `loadfid` does, so `loadfid(savefid(...))` round-trips.

### Peaks and assignments

- `findpeak` refines a single peak near a given position, and `findnearest` finds the closest peak
  in a list.
- `register!` and `register` translate a peak list onto a spectrum by a global search.
- `savepeaks` and `loadpeaks` read and write NMRAnalysis.jl `peaklist.csv` files and Sparky lists.
  A peak list is also a Tables.jl table, so `DataFrame(pks)` works.
- `loadbmrb` reads assigned chemical shifts from an NMR-STAR file or by BMRB ID, and `peaklist`
  builds peak lists from them by atom name (with residue offsets) or by nucleus.
- Peak lists can be plotted with Plots (as lines in 1D, crosses in 2D) and with Makie (`nmrplot`
  and `nmrplot!`), with or without a spectrum.

### Other changes

- Plots: the normalisation logic of the contour and pseudo-2D recipes has been consolidated. In
  the recipe for several pseudo-2D datasets, `normalize=refspec` previously divided every dataset
  by the reference's scale, leaving them unnormalised relative to each other; each is now divided by
  its own scale, as in the other recipes. An invalid `normalize` value now throws an
  `ArgumentError` in the 1D recipes too.
- Data loading is included in the precompilation workload, reducing the time to first load.
