# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Common Development Commands

**Testing:**
- `julia -e 'using Pkg; Pkg.test()'` — run all tests using SafeTestsets
- Test files: `test/nmrbase_test.jl`, `test/nmrio_test.jl`, `test/plotsext_test.jl`

**Documentation:**
- `julia make-local-docs.jl` — build documentation locally (installs package in dev mode and runs docs/make.jl)
- Documentation uses Documenter.jl with visual regression tests for plots
- Do not add "Reference" sections (with `@docs` or `@autodocs` blocks) to doc pages — all API reference lives in `docs/src/api.md` only
- Documenter.jl resolves `@ref` links in the context of the module where a docstring is defined. When a docstring in `NMRBase` references a function in `NMRIO`, the unqualified `@ref` fails. Use the fully-qualified form: `` [`name`](@ref NMRTools.NMRIO.name) ``

**Package Management:**
- Standard Julia package — use `] dev .` to develop locally
- Extensions system: PlotsExt.jl provides plotting functionality when Plots.jl is loaded

## Architecture Overview

NMRTools.jl is a Julia library for NMR spectroscopy data. It encodes the structure of NMR data in the type system: multicomplex algebra for quadrature detection, named dimensions for physical axes, and metadata for acquisition parameters.

### Module structure

- **NMRTools.jl** (main): re-exports NMRBase, NMRIO, and MulticomplexNumbers
- **NMRBase/**: core data structures, dimension types, metadata, window functions, lineshapes, peaks, referencing
- **NMRIO/**: file I/O (Bruker, nmrPipe, UCSF, JDX), annotation parsing, sample metadata
- **ext/PlotsExt.jl**: plotting recipes loaded when Plots.jl is available

### Key types

- **`NMRData <: AbstractDimArray`**: main data container. Wraps arrays with named dimensions and metadata. Built on DimensionalData.jl.
- **`Multicomplex{T,N,C}`**: element type for NMR data. From MulticomplexNumbers.jl. A 2D experiment has `Multicomplex{Float64,2,4}` elements (bicomplex, 4 components). `im1` is the conventional complex unit — never use Julia's `im`.
- **Dimension hierarchy**: `FrequencyDimension` (F1Dim–F4Dim), `TimeDimension` (T1Dim–T4Dim, TrelaxDim, TkinDim), `GradientDimension` (G1Dim–G4Dim), `UnknownDimension` (X1Dim–X4Dim), `OffsetDimension` (OffsetDim), `FieldDimension` (SpinlockDim)
- **`WindowFunction`** subtypes: `ExponentialWindow`, `CosWindow`, `Cos2Window`, `SineWindow`, `GaussWindow`, `LorentzToGaussWindow`, `NullWindow`, `UnknownWindow`

### Critical metadata: `:mcindex`

Every `TimeDimension` and `FrequencyDimension` has an `:mcindex` metadata entry that links the dimension to its imaginary unit in the multicomplex representation. This is set once by `loadfid` and never modified:

| Dimension | `:mcindex` | Imaginary unit |
|---|---|---|
| T1Dim / F1Dim | `1` | `im1` |
| T2Dim / F2Dim | `2` | `im2` |
| T3Dim / F3Dim | `3` | `im3` |
| T4Dim / F4Dim | `4` | `im4` |
| TrelaxDim, X2Dim, etc. | `nothing` | none (not Fourier transformed) |

Processing functions (`ft`, `phase`, `apodize`) read `:mcindex` from the target dimension to determine which imaginary unit to operate on. The `fft!`/`ifft!` functions in MulticomplexNumbers.jl accept the mcindex directly as their `unit` argument.

### Data representation

The element type of `NMRData` is `Multicomplex{Float64,N,C}` throughout the entire processing pipeline — from `loadfid` through apodization, zero-filling, Fourier transform, and phase correction. The element type never changes during processing. Only `real(spec)` at the very end extracts the real part and converts to `Float64`.

Even 1D data uses `Multicomplex{Float64,1,2}` (first-order, isomorphic to Complex but using `im1` not `im`). This keeps the type system uniform. Use `complex(spec)` for explicit conversion to Julia `Complex` when needed.

### File format support

- **Bruker processed** (pdata): `loadnmr("experiment/1/pdata/1")` or just `loadnmr("experiment/1")` — 1D, 2D, 3D
- **nmrPipe**: `loadnmr("test.ft2")` — 1D, 2D, 3D (including plane series)
- **UCSF/Sparky**: `loadnmr("spectrum.ucsf")`
- **`loadnmr`** auto-detects format from filename/path

Reading time-domain data (under development):
- **Bruker raw** (fid/ser): `loadfid("experiment/1")` — under development
- Digital filter: `loadfid` keeps the `GRPDLY` build-up points on a time axis starting at `-GRPDLY` dwells and records `:grpdly`; `ft` removes the delay with a linear phase ramp after the FFT (as TopSpin/NMRPipe — matches TopSpin 1r to ~1e-8). Dropping the points instead gives the edge "smile". `apodize` leaves t < 0 unwindowed; `linearpredict` works from t = 0. Baseopt (`DIGMOD = 3`): t = 0 point not scaled. `DSPFVS` outside 10–23 (e.g. 0 after TopSpin's `cryoproc1d`, which rewrites the FID with its start removed): no delay, warn.

Writing Bruker data:
- **`savenmr(path, spec; procno=1)`**: frequency-domain data → `pdata/<procno>/1r`, `2rr`, … (int32 + NC_proc)
- **`savefid(path, fid)`**: time-domain data → raw `fid`/`ser`, undoing everything `loadfid` did (first-point scaling, group delay, quadrature recombination, axis order); `loadfid(savefid(...))` round-trips
- Both copy parameter files from the source experiment (`:experimentfolder` metadata) or generate minimal ones

Sign convention: NMRTools spectra are the complex conjugate of TopSpin's in every imaginary unit (real parts agree). `phase` uses TopSpin's sign, so `PHC0`/`PHC1` values carry over directly. `loadnmr(...; allcomponents=true)` and `savenmr` conjugate imaginary components when reading/writing pdata.

### Annotation system (under development)

Pulse programme annotations (`;@ ` prefixed YAML in pulse programmes) are parsed automatically by `annotate()` during loading. They resolve parameter references against acqus, generate programmatic lists, and apply semantic dimension types. Schema v0.0.2. See `ANNOTATIONS.md` for the schema.

### Processing pipeline (under development)

Target API uses `|>` piping with composable functions:

```julia
loadfid("1/ser") |> apodize |> zerofill |> ft |>
    phase(F1Dim, -50) |> phase(F2Dim, 0) |> real
```

Processing functions dispatch on dimension types — `TimeDimension`s are processed, others are skipped. All functions accept optional dimension arguments for per-dimension control.

`autoprocess(fid; window, zerofill, solvent, phase, rebuild)` chains solvent filter (¹H in H₂O), apodize (`Cos2Window` default), zerofill, ft and autophase for 1D/2D; per-dimension settings via `Dict`. `loadfid` scales raw data by `2^NC` so intensities match TopSpin. `phase`, `extract` and `ft` refresh `:noise`. `Cos2Window` (alias `Cos²Window`).

Automatic phasing: `estimatepoles` (band-wise HSVD → `Vector{Pole}`), `estimatephase` (doubled-angle circular regression of multiplet-grouped pole phases), `autophase` (estimate + apply). `phase` accumulates the total correction in the dimension's `:phase` metadata.

## Style Conventions

Follow Julia Base style throughout:

- **Public function names**: lowercase, no underscores (`loadnmr`, `apodize`, `zerofill`, `estimatenoise`, `detectpeaks`, `lineshape`). Avoid underscores unless necessary for readability of long compound names.
- **Internal functions**: underscore prefix (`_buildfreqdim`, `_readfid`, `_resolve_recursive!`)
- **Types**: CamelCase (`NMRData`, `FrequencyDimension`, `ExponentialWindow`)
- **Constants**: UPPER_SNAKE_CASE for true constants, CamelCase for type-like constants (`im1`, `im2` follow MulticomplexNumbers convention)

Prefer multiple dispatch over conditional logic. Use positional arguments for the primary operation, keyword arguments for options. Argument ordering follows Julia Base convention: function target first, then modifiers.

## Error Handling

- Use Julia standard exceptions where they apply: `ArgumentError` for bad arguments, `BoundsError` for out-of-range indices, `DimensionMismatch` for size mismatches.
- Use `NMRToolsError` only for domain-specific errors: unknown file format, missing `:mcindex`, unsupported NMR operation.
- Never use bare `error()` — always throw a typed exception.
- Never silently swallow exceptions in `catch` blocks — at minimum `@debug` the error.

## Logging

- **`@debug`**: internal tracing, parsing details, metadata construction. Annotation parsing failures go here.
- **`@info`**: significant user-facing operations (referencing applied, data loaded, format detected).
- **`@warn`**: things the user should know about and may need to act on (unsupported schema version, missing calibration data, data size mismatches, silent assumptions like default gmax).
- Avoid `@warn` for best-effort parsing that commonly fails gracefully — use `@debug` instead.

## Multicomplex Conventions

- `im1` is the conventional NMR complex unit. Never use Julia's `im`.
- A `Multicomplex{Float64,N,C}` has `N` imaginary units (`im1` through `imN`) and `C = 2^N` real components.
- Component ordering: for bicomplex (N=2), components are `(rr, ir, ri, ii)` representing `rr + ir*im1 + ri*im2 + ii*im1*im2`.
- `real(m)` strips the highest-order imaginary unit. `realest(m)` extracts the purely real component.
- `fft!(A, unit, dims)` Fourier transforms array `A` along array dimensions `dims`, treating imaginary unit `unit` as the complex `i` for the transform.

## Key Dependencies

- **MulticomplexNumbers.jl** (`waudbylab/MulticomplexNumbers.jl`): multicomplex number types with FFTW extension for `fft!`/`ifft!`
- **DimensionalData.jl**: named dimension arrays, selectors, metadata infrastructure
- **Plots.jl** (weak dependency): plotting recipes via extension
- **FFTW.jl** (via MulticomplexNumbers): Fourier transforms

## Testing Patterns

- SafeTestsets for isolated test environments
- Visual regression tests for plotting (VisualRegressionTests.jl with reference PNGs)
- Test data: real NMR datasets bundled as Julia Artifacts (Bruker and nmrPipe formats)
- Example data accessed via `exampledata("2D_HN")` etc.
- The precompile workload (`src/precompile.jl`) loads `test/test-data/samples/2`, `test/test-data/19f-cest-ts3` and `test/test-data/15n-trosy-ts3`, since the example artifacts are lazy and would need a download. Keep these datasets in place, or update the workload if they move.