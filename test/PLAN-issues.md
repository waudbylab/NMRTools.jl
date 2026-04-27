# NMRTools.jl GitHub Issues

Issues grouped by milestone, roughly in dependency order within each group.


## Milestone: Test Coverage

**8. Tests for sumexpts**
File I/O, weighted summation, metadata updates.

**9. Tests for remaining window functions**
All `apod` methods: `GaussWindow`, `GeneralGaussWindow`, `LorentzToGaussWindow`, `GeneralSineWindow`.

**10. Edge case tests for format readers**
Malformed acqus, missing parameters, unusual dimensionalities.


## Milestone: MulticomplexNumbers.jl

**17. Add ifft! to FFTW extension**
`ifft!(A, unit)` and `ifft!(A, unit, dims)`. Mirror `fft!` pattern. Required for inverse FT and IST.


## Milestone: Raw Time-Domain Loading

**18. Implement loadfid**
Extend `getformat` for fid/ser. Return `Multicomplex` element type. Set `:mcindex` per dimension. Automatic quadrature recombination (States/States-TPPI/echo-anti-echo from FnMODE or annotations). Automatic digital filter removal (opt-out). First-point scaling (default 0.5, overridable). NUS detection. Call `annotate`/`addsampleinfo`.

**19. Digital filter lookup table**
DECIM/DSPFVS to group delay for legacy spectrometers. Write from scratch citing Westler and Abildgaard.

**20. Tests for loadfid**
1D fid, 2D ser with each quadrature scheme, pseudo-2D, digital filter vs TopSpin reference, NUS detection, round-trip vs `loadpdata`.


## Milestone: Processing Pipeline

**21. Implement apodize**
`apodize(spec, [dim]; window=:default, inverse=false)`. Operates on `TimeDimension`s. Curried form. Wraps existing `WindowFunction`/`apod`.

**22. Implement zerofill**
`zerofill(spec, [dim]; factor=2)`. Operates on `TimeDimension`s. Curried form.

**23. Implement ft**
`ft(spec, [dim]; inverse=false)`. Reads `:mcindex`, calls `fft!`/`ifft!`, applies fftshift, builds new dimension with chemical shift axis. No-argument form processes all time dimensions in order.

**24. Implement phase**
`phase(spec, dim, p0, p1=0; pivot=nothing)`. Reads `:mcindex` for imaginary unit. Works in time or frequency domain. Curried form.

**25. Implement extract**
`extract(dim, range)` returns closure for piping. Support multiple dimensions via pairs.

**26. Extend Base.real, abs, imag, complex for NMRData**
`real` applies `realest`, rebuilds as `Float64`. `complex` converts first-order multicomplex.

**27. Collection mapping for processing functions**
Accept `AbstractDict{Symbol, <:NMRData}`, map over values.

**28. Tests for processing pipeline**
Round-trip vs `loadpdata`. Phase 360 identity. FT inverse identity. Extract equivalence. Type conversions.


## Milestone: Write Support

**29. Implement savenmr with format inference**
`savenmr(filename, spec; format=nothing)`. Infer from extension: `.ft1`/`.ft2`/`.ft3` for nmrPipe, `.nmr`/`.jld2` for JLD2, pdata path for Bruker. Optional `format` keyword override.

**30. nmrPipe write**
Populate header from metadata, write Float32. Support 1D, 2D, 3D.

**31. JLD2 serialisation**
Save/load NMRData. Schema version for compatibility. Extend `getformat`.

**32. Bruker pdata write (lower priority)**
Use existing `reordersubmatrix(reverse=true)`. Generate procs from metadata.

**33. Tests for write support**
nmrPipe round-trip. JLD2 round-trip. Bruker pdata round-trip if implemented.
