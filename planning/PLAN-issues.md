# NMRTools.jl outstanding issues

Work remaining after v0.3.0, grouped by milestone and roughly in dependency order within each group.
The completed milestones (MulticomplexNumbers `ifft!`, raw time-domain loading, the processing
pipeline and Bruker write support) have been removed. See `PLAN-development.md` for the designs.


## Milestone: Test coverage

**Tests for sumexpts**
File I/O, weighted summation, metadata updates.

**Tests for remaining window functions**
All `apod` methods: `GaussWindow`, `GeneralGaussWindow`, `LorentzToGaussWindow`, `GeneralSineWindow`.

**Edge case tests for format readers**
Malformed acqus, missing parameters, unusual dimensionalities.

**Test loadbmrb against the live BMRB** (GitHub #78)


## Milestone: Raw time-domain loading

**NUS detection and schedule loading**
Detect NUS from `FnTYPE` or a `nuslist` file, place acquired FIDs in a zero-filled grid, record
`:nus` and `:nuslist` metadata.

**4D data**
Extend `loadfid`, `savefid` and `savenmr` beyond 3D.

**Legacy digital filter lookup table**
`(DSPFVS, DECIM)` to group delay for data without `GRPDLY`. Write from scratch citing Westler and
Abildgaard.


## Milestone: Write support

**Format inference in savenmr**
`savenmr(path, spec; format=nothing)` infers the format from the extension: `.ft1`/`.ft2`/`.ft3` for
nmrPipe, `.nmr`/`.jld2` for JLD2, otherwise Bruker pdata as now.

**nmrPipe write**
Populate header from metadata, write Float32. Support 1D, 2D, 3D.

**JLD2 serialisation**
Save/load NMRData. Schema version for compatibility. Extend `getformat`.

**Tests for write support**
nmrPipe round trip. JLD2 round trip.


## Milestone: NUS reconstruction

**IST**
`ist(spec, dims...; iterations, threshold)` with a curried form, replacing
`apodize |> zerofill |> ft` for NUS dimensions.


## Milestone: Plots extension

**Pseudo-2D recipe defaults**
Reconsider the `:path3d` waterfall default; remove or fix the unused `usegradient` attribute.


## Milestone: Code hygiene

**Silent catch in NMRBase/metadata.jl**
Add at least a `@debug` with the exception, as the error-handling conventions require.
