# NMRTools.jl Development Plan

## Design Principles

1. **Domain modelling first.** Types encode the structure of NMR data: multicomplex algebra for quadrature, named dimensions for physical axes, metadata for acquisition parameters. Processing functions dispatch on these types.

2. **Julian idioms.** Extend `Base` and standard library functions where appropriate (`real`, `complex`, `abs`, `fft!`). Public function names are lowercase without underscores, following Julia Base convention (`loadfid`, `apodize`, `zerofill`, `lineshape`). Types are CamelCase. Underscore-prefixed names are internal only. Use Julia standard exception types where applicable; reserve custom exceptions for domain-specific errors.

3. **Multicomplex numbers are the data representation.** `im1` is the conventional complex unit. Element types never change during processing — only `real(spec)` at the end of a pipeline strips imaginary components. There is no mixing with Julia's built-in `Complex`; use `complex(spec)` to convert a first-order multicomplex to `Complex` when needed.

4. **`:mcindex` metadata links dimensions to imaginary units.** Set once by `loadfid`, never modified, carried through all operations. Every processing function reads `:mcindex` from the target dimension's metadata to determine which imaginary unit to operate on.

5. **Composable pipeline.** Processing functions take and return `NMRData`. Curried forms enable `|>` piping. No-argument forms apply sensible defaults to all applicable dimensions. Explicit dimension arguments available when control is needed.

---

## 1. Repository Hygiene (~1 day)

- Add `CONTRIBUTING.md` and `CODE_OF_CONDUCT.md`.
- Ensure GitHub Actions CI with Codecov badge in README.
- Verify `CITATION.cff` is complete.
- Remove "under active development" warning from docs homepage.
- Audit README: install instructions, minimal example, doc links.

---

## 2. Error Handling and Logging Audit (~1 day)

### 2.1 Exception types

Currently everything throws `NMRToolsError`. Adopt Julian convention: use standard exception types where they apply, reserve `NMRToolsError` for genuinely domain-specific errors.

| Situation | Current | Should be |
|---|---|---|
| Dimension index out of range | `NMRToolsError` | `BoundsError` |
| Wrong dimension type for operation | `NMRToolsError` | `ArgumentError` |
| Incompatible data sizes | `NMRToolsError` | `DimensionMismatch` |
| Bad keyword value | `error("...")` | `ArgumentError` |
| Unknown NMR file format | `NMRToolsError` | `NMRToolsError` (domain-specific, keep) |
| Cannot FT: no `:mcindex` | `NMRToolsError` | `NMRToolsError` (domain-specific, keep) |
| File not found / can't load | `NMRToolsError` | `SystemError` or `ArgumentError` |

Replace all bare `error("...")` calls (in `nmrdata.jl:61`, Plots extension) with typed exceptions.

Fix the silent `catch e` in `loadnmr.jl:73` — at minimum `@debug "Could not compute field strength" exception=e`.

### 2.2 Logging levels

Audit all `@warn`, `@info`, `@debug` usage for consistency:

**`@debug`** — internal tracing, parsing details. Currently correct usage in `dimensions.jl`, `acqus.jl`. Demote most `@warn` in annotation parsing (`annotation.jl`) to `@debug`: these are best-effort parsing messages that are noisy during normal use (missing counter values, missing programmatic list parameters, failed YAML parse).

**`@info`** — significant user-facing operations. Currently correct in `reference.jl` (reporting referencing changes). Use for `loadfid` summary (format detected, dimensions identified, digital filter removed).

**`@warn`** — things the user should know about and may need to act on. Keep for: unsupported schema version, missing gmax assumption, data/annotation size mismatches, inability to locate acqus file, silent data truncation, experiment stack inconsistencies (`ns`/`rg` mismatch in `sumexpts`).

**`@error`** — reserved for cases where something has gone wrong but execution can continue (rare in NMRTools; prefer throwing exceptions).

---

## 3. Test Coverage (ongoing)

Priority gaps:

- **Annotation system** (`annotation.jl`): parsing, parameter resolution, programmatic lists, dimension application.
- **Multicomplex loading** (`loadpdata` with `allcomponents=true`): component assembly, correct `Multicomplex` construction.
- **`sumexpts`**: file I/O and weighted summation.
- **Window functions**: all `apod` methods including `GaussWindow`, `GeneralSineWindow`.
- **Edge cases**: malformed acqus, missing parameters, unusual dimensionalities.

Target: >80% line coverage, with particular depth on multicomplex and annotation paths.

---

## 4. Plots Extension Fixes (~2–3 days)

### 4.1 Contour colour handling

The current logic extracts hue via `sequential_palette(hue(convert(HSV, basecolor)))`, which destroys saturation and lightness and produces red for any achromatic colour (black, white, grey).

Replace with explicit `poscolor`/`negcolor` keyword arguments:

```julia
# New keywords
poscolor   # colour for positive contours (default: derived from seriescolor)
negcolor   # colour for negative contours (default: lighter version of poscolor)
negcontours  # Bool, whether to show negative contours at all (default: true)
```

Automatic derivation works in HSVA space preserving the original colour's properties:

```julia
hsv = convert(HSVA, poscolor)
negcolor = HSVA(hsv.h, hsv.s * 0.5, min(1.0, hsv.v + 0.3), hsv.alpha)
```

This handles black correctly (gives dark grey), preserves saturated colours, and respects user choices. The `seriescolor` / `c` keyword continues to work as the base for automatic derivation.

For overlay plots (vector of 2D spectra), accept vectors: `poscolor=[:blue, :red, :black]`.

Add `negcontours=false` option to suppress negative contours entirely.

### 4.2 Normalisation logic

The single-spectrum 2D recipe has `normalize == true` and `normalize == false` branches that compute the same thing (`σ = dfwd[:noise]`). Fix: `normalize == false` should set `σ` to some absolute value or skip normalisation-based contour levels entirely.

Deduplicate the normalisation logic between single-spectrum and vector recipes — extract into a helper function.

### 4.3 Pseudo-2D recipe

The default `seriestype --> :path3d` waterfall is unusual. Consider changing the default to a heatmap or stacked 1D overlay, which are more standard for NMR pseudo-2D data. Keep 3D waterfall available via explicit `seriestype=:path3d`.

The `usegradient` attribute is read but never works (deleted from attributes before the dispatch path reaches it). Either fix the dispatch or remove the dead code.

### 4.4 Broken recipes

The vector-of-pseudo-2D recipe references `HasPseudoDimension` which doesn't exist — this will error at runtime. Either implement properly or replace with a clear error message.

### 4.5 Multicomplex data support

Add methods for plotting `NMRData` with `Multicomplex` element types. At minimum, take `realest` before delegating to existing recipes:

```julia
@recipe function f(A::NMRData{<:Multicomplex, N}) where N
    realest.(parent(A))  # or rebuild with realest data
end
```

This prevents plotting from breaking once `loadfid` returns multicomplex data.

### 4.6 Error handling in Plots extension

Replace bare `error("normalize must be true, false or a reference spectrum")` with `throw(ArgumentError(...))`.

### 4.7 Tests

- Test contour colours with black, white, grey, and saturated inputs.
- Test `poscolor`/`negcolor` explicit specification.
- Test `negcontours=false`.
- Update visual regression tests as needed.

---

## 5. Raw Time-Domain Loading (~2–3 weeks)

### 5.1 `loadfid`

New entry point for raw Bruker fid/ser files. Extend `getformat` to recognise `fid` and `ser` files.

**Signature:**

```julia
loadfid(filename; removefilter=true, firstpoint=nothing)
```

**Data representation by experiment type:**

| Experiment | Element type | Dimensions |
|---|---|---|
| 1D | `Multicomplex{Float64,1,2}` | `(T1Dim,)` |
| Pseudo-2D | `Multicomplex{Float64,1,2}` | `(T1Dim, X2Dim)` |
| 2D | `Multicomplex{Float64,2,4}` | `(T1Dim, T2Dim)` |
| Pseudo-3D | `Multicomplex{Float64,2,4}` | `(T1Dim, T2Dim, X3Dim)` |
| 3D | `Multicomplex{Float64,3,8}` | `(T1Dim, T2Dim, T3Dim)` |
| 4D | `Multicomplex{Float64,4,16}` | `(T1Dim, T2Dim, T3Dim, T4Dim)` |

Note: even 1D data uses `Multicomplex{Float64,1,2}` (first-order multicomplex, isomorphic to `Complex`), not Julia `Complex`. This keeps the type system uniform — `im1` is always the imaginary unit.

**`:mcindex` assignment:**

| Dimension | `:mcindex` | Used by `fft!` as `unit` argument |
|---|---|---|
| T1Dim (direct) | `1` | `fft!(data, 1, arraydim)` |
| T2Dim (1st indirect) | `2` | `fft!(data, 2, arraydim)` |
| T3Dim (2nd indirect) | `3` | `fft!(data, 3, arraydim)` |
| T4Dim (3rd indirect) | `4` | `fft!(data, 4, arraydim)` |
| TrelaxDim, X2Dim, etc. | `nothing` | not FT'd |

**TimeDimension metadata set by `loadfid`:**

```
:mcindex     => 2               # imaginary unit index
:label       => "15N"           # nucleus label
:nucleus     => N15             # Nucleus enum
:bf          => 60817738.0      # base frequency (Hz)
:swhz        => 2500.0          # spectral width (Hz)
:swppm       => 41.1            # spectral width (ppm)
:offsetppm   => 118.0           # carrier offset (ppm)
:offsethz    => ...             # carrier offset (Hz)
:sf          => ...             # carrier frequency (Hz)
:td          => 128             # complex points acquired
:window      => nothing         # no window applied yet
:quadrature  => :states_tppi    # quadrature scheme used
```

**Quadrature recombination** is handled inside `loadfid`. Read `FnMODE` from acqus (or annotations). Apply the appropriate recombination (States, States-TPPI, echo-anti-echo) including sign corrections (negate imaginary, alternating sign). Store scheme in metadata for reference. The user never sees interleaved FIDs.

**Digital filter removal** is automatic by default (`removefilter=true`). Modern path: read `GRPDLY` from acqus, apply circular shift + first-order phase. Legacy path: lookup table from `(DSPFVS, DECIM)` → group delay, written from scratch citing Westler & Abildgaard. Override with `removefilter=false`.

**First-point scaling** defaults to 0.5 for all time dimensions. Override per dimension: `firstpoint=Dict(T2Dim => 1.0)`. When annotations specify first-order phase corrections, infer the correct scaling automatically.

**NUS detection:** if `FnTYPE` indicates NUS or `nuslist` file present, load schedule, place acquired FIDs at correct positions in zero-filled grid, store schedule in metadata:

```
T2Dim metadata[:nus] = true
T2Dim metadata[:nuslist] = [3, 7, 12, ...]
```

**Annotations and sample metadata:** call `annotate(spec)` and `addsampleinfo(spec)` as `loadnmr` already does. Annotations set semantic dimension types automatically.

### 5.2 Binary reader

Implement `_readfid(filename, acqusdict)` (internal). Handle endianness (`BYTORDA`), data type (`DTYPA`), direct-dimension padding, interleaved real/imaginary → multicomplex assembly.

### 5.3 Tests

- 1D fid loading, known values.
- 2D ser with States, States-TPPI, echo-anti-echo.
- Pseudo-2D (relaxation series).
- Digital filter removal vs TopSpin reference.
- NUS detection and schedule loading.
- Round-trip: `loadfid` → process → compare with `loadpdata`.

---

## 6. Processing Pipeline (~2–3 weeks)

### 6.1 Core functions

All processing functions:
- Accept `NMRData` and optional dimension identifier (dimension type, integer, or `Nucleus` — error if ambiguous).
- No-argument dimension defaults: all `TimeDimension`s for `apodize`/`zerofill`/`ft`.
- Return new `NMRData` (non-mutating). Provide `!` variants where performance matters.
- Curried forms for `|>` piping.
- Element type stays `Multicomplex` throughout. Only `real` converts.

**`apodize(spec, [dim]; window=:default, inverse=false)`**

Apply window function along specified `TimeDimension`(s). `:default` reads from dimension metadata (`:window`). `inverse=true` divides by the window function. The existing `WindowFunction` type hierarchy and `apod()` functions provide the foundation.

**`zerofill(spec, [dim]; factor=2)`**

Zero-fill `TimeDimension`(s) by specified factor. Updates dimension length and `:tdzf` metadata.

**`ft(spec, [dim]; inverse=false)`**

Fourier transform along specified dimension(s).

Implementation for a single dimension `di`:

```julia
function ft(spec::NMRData, dim; inverse=false)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    mcindex = metadata(d, :mcindex)

    isnothing(mcindex) && throw(NMRToolsError(
        "Cannot FT $(typeof(d)): no :mcindex"))

    newdata = copy(parent(spec))

    if inverse
        newdata = ifftshift(newdata, di)
        ifft!(newdata, mcindex, di)
        newdim = _buildtimedim(d, di, size(newdata, di))
    else
        fft!(newdata, mcindex, di)
        newdata = fftshift(newdata, di)
        newdim = _buildfreqdim(d, di, size(newdata, di))
    end

    newdims = _replacedimtuple(dims(spec), di, newdim)
    rebuild(spec, newdata, newdims)
end
```

Where `_buildfreqdim` constructs the new `FrequencyDimension` (F1Dim/F2Dim/F3Dim/F4Dim) with:
- Chemical shift axis computed from `:bf`, `:swhz`, `:offsetppm`
- All metadata carried over from the time dimension, including `:mcindex`
- `:units => "ppm"`, `:tdzf` and `:npoints` set to new size, `:region => missing`

And `_buildtimedim` does the reverse for inverse FT.

No-argument `ft(spec)` processes all `TimeDimension`s in dimension order (T1Dim first = direct first).

`fft!` and `ifft!` come from MulticomplexNumbers.jl's FFTW extension. The `unit` argument is read directly from `:mcindex`. The `dims` argument is the array dimension index `di`.

**`phase(spec, dim, p0, p1=0; pivot=nothing)`**

Phase correction along a specified dimension. Reads `:mcindex` from the dimension to determine which imaginary unit carries the phase. Works identically in time or frequency domain — multiplication by a multicomplex phase factor `exp(im_k * θ)` where `k` is the mcindex.

- `p0` in degrees (positional).
- `p1` in degrees (positional, default 0).
- `pivot` in ppm (frequency domain) or seconds (time domain). Default: carrier frequency.

Curried: `phase(dim, p0, p1=0; pivot=nothing)` returns a closure.

**`extract(dim, range)`**

Pipeline-compatible region selection:

```julia
extract(dim, range) = spec -> spec[dim(range)]
```

Multiple dimensions: `extract(F1Dim => 6.0..10.0, F2Dim => 110.0..130.0)`.

**`real(spec::NMRData)`**

Extend `Base.real`. Applies `realest` elementwise and rebuilds with `Float64` element type:

```julia
Base.real(spec::NMRData{<:Multicomplex}) = rebuild(spec, realest.(parent(spec)))
```

Similarly extend `Base.abs` and `Base.imag`.

**`Base.complex(spec::NMRData{<:Multicomplex{T,1,2}})`** — convert first-order multicomplex to Julia `Complex{T}`.

### 6.2 Collection mapping

All processing functions accept `AbstractDict{Symbol, <:NMRData}` and map over values:

```julia
apodize(d::AbstractDict{Symbol, <:NMRData}; kw...) =
    Dict(k => apodize(v; kw...) for (k, v) in d)
```

### 6.3 Target syntax

```julia
# 1D
loadfid("1/fid") |> apodize |> zerofill |> ft |> phase(F1Dim, -50) |> real

# 2D HSQC
loadfid("1/ser") |> apodize |> zerofill |> ft |>
    phase(F1Dim, -50) |> phase(F2Dim, 0) |> real

# 2D with per-dimension control and extraction
loadfid("1/ser") |>
    apodize(T1Dim; window=ExponentialWindow(5.0)) |>
    apodize(T2Dim; window=CosWindow()) |>
    zerofill |> ft(T1Dim) |>
    extract(F1Dim, 6.0 .. 10.0) |>
    ft(T2Dim) |>
    phase(F1Dim, -50) |> phase(F2Dim, 0) |> real

# Pseudo-2D relaxation series
loadfid("1/ser") |> apodize |> zerofill |> ft |>
    phase(F1Dim, -50) |> real

# 3D HNCO with progressive extraction
loadfid("1/ser") |>
    apodize(T1Dim; window=ExponentialWindow(5.0)) |>
    zerofill(T1Dim) |> ft(T1Dim) |>
    extract(F1Dim, 6.0 .. 10.0) |>
    apodize(T2Dim; window=CosWindow()) |>
    zerofill(T2Dim) |> ft(T2Dim) |>
    extract(F2Dim, 110.0 .. 130.0) |>
    apodize(T3Dim; window=Cos²Window()) |>
    zerofill(T3Dim) |> ft(T3Dim) |>
    extract(F3Dim, 170.0 .. 180.0) |>
    phase(F1Dim, -50) |> phase(F2Dim, 0) |> phase(F3Dim, 0) |> real
```

### 6.4 Tests

- Round-trip: `loadfid` → process → compare with `loadpdata`.
- `phase(dim, 360) ≈ identity`.
- Multicomplex FT: verify against reference.
- `extract` within pipeline vs post-hoc slicing.
- `real`/`abs`/`imag`: type conversion, data correctness.
- `ft(inverse=true)` ∘ `ft` ≈ identity (up to normalisation).

---

## 7. Write Support (~1 week)

### 7.1 Unified `savenmr`

Single save function that infers format from file extension:

```julia
savenmr("output.ft2", spec)      # nmrPipe format (from .ft1/.ft2/.ft3 extension)
savenmr("output.nmr", spec)      # JLD2 serialisation (from .nmr/.jld2 extension)
savenmr("output/pdata/1", spec)  # Bruker pdata (from path structure)
savenmr("output.ucsf", spec)     # UCSF/Sparky (future)
```

Optional keyword to override format detection: `savenmr("output.dat", spec; format=:nmrpipe)`.

Corresponding loading already handled by `loadnmr` with `getformat` dispatcher — extend `getformat` to recognise `.nmr`/`.jld2` for serialised data.

### 7.2 nmrPipe write

Populate header from NMRData metadata (the header specification is already fully parsed in `parsenmrpipeheader`). Write header + Float32 data. Support 1D, 2D (single file), 3D (plane series).

### 7.3 JLD2 serialisation

Workflow persistence for Julia-to-Julia use. Store schema version for forward compatibility. Not interoperable with other software.

### 7.4 Bruker pdata write (lower priority)

Submatrix reordering already exists (`reordersubmatrix` with `reverse=true`). Need procs file generation from metadata.

### 7.5 Tests

- nmrPipe round-trip: write, read back, verify.
- JLD2 round-trip with various NMRData types.
- Bruker pdata round-trip if implemented.

---

## 8. NUS Reconstruction (~1–2 weeks, post Phase 6)

### 8.1 IST

```julia
ist(spec, dim; iterations=200, threshold=0.98)
ist(spec, dim1, dim2; iterations=300)    # joint reconstruction
```

Curried: `ist(dim; kw...)` for piping.

IST loop uses the same `fft!`/`ifft!` and `apodize` internals as the main pipeline. Reads `:mcindex` from dimensions, `:nuslist` and `:window` from metadata.

`ist` replaces `apodize |> zerofill |> ft` for NUS dimensions. Outputs frequency-domain data with `FrequencyDimension` types.

Direct dimension can be extracted before IST. NUS dimensions must remain full-width during reconstruction.

**Workflow:**

```julia
# 2D NUS HSQC
loadfid("1/ser") |>
    apodize(T1Dim) |> zerofill(T1Dim) |> ft(T1Dim) |>
    extract(F1Dim, 6.0 .. 10.0) |>
    ist(T2Dim; iterations=200) |>
    phase(F1Dim, -50) |> phase(F2Dim, 0) |> real

# 3D NUS HNCO — joint reconstruction of both indirect dimensions
loadfid("1/ser") |>
    apodize(T1Dim) |> zerofill(T1Dim) |> ft(T1Dim) |>
    extract(F1Dim, 6.0 .. 10.0) |>
    ist(T2Dim, T3Dim; iterations=300) |>
    phase(F1Dim, -50) |> phase(F2Dim, 0) |> phase(F3Dim, 0) |> real
```

### 8.2 Tests

- Retrospectively undersampled 2D, compare with fully sampled.
- 100% sampling reproduces standard FT.
- Schedule loading and mask generation.

---

## 9. Annotation Schema Extensions (parallel)

### 9.1 Quadrature annotations

```yaml
;@ dimensions:
;@   - f1
;@   - name: f2
;@     quadrature:
;@       scheme: echo-anti-echo
;@       negate_imaginary: true
;@       alternate_sign: false
```

Read by `loadfid` for automatic recombination. Falls back to acqus `FnMODE`.

### 9.2 Processing hints (forward-looking)

```yaml
;@ processing:
;@   - dimension: f1
;@     window: {type: exponential, lb: 5.0}
;@     phase: {p0: 0, p1: 0}
;@   - dimension: f2
;@     window: {type: cosine_squared}
;@     phase: {p0: 0, p1: 180}
;@     firstpoint: 0.5
```

Design the schema now. Implement `autoprocess` later.

### 9.3 Phase cycle annotations (forward-looking)

```yaml
;@ phasecycle:
;@   steps: 4
;@   receiver: [0, 180, 0, 180]
;@   experiments:
;@     zqc: {weights: [1, -1, 1, -1]}
;@     dqc: {weights: [1, 1, -1, -1]}
```

Design the schema now. Implement `loadfid(; split=true)` later.

---

## 10. JOSS Paper

### Minimum viable for submission

- `loadfid` with multicomplex representation (Phase 5).
- Processing pipeline for 1D and 2D (Phase 6).
- `savenmr` with nmrPipe format (Phase 7).
- Test coverage >80% (Phase 3).
- Repository hygiene (Phase 1).
- Error handling / logging audit (Phase 2).
- Plots extension fixes (Phase 4).
- Updated documentation with processing tutorials.

### Paper structure

- **Summary**: Julia library encoding the algebraic and dimensional structure of NMR data in the type system.
- **Statement of need**: NMR data has three layers of structure (multicomplex quadrature, named physical dimensions, acquisition metadata) that existing tools represent only by convention. NMRTools.jl makes them explicit.
- **State of the field**: nmrglue (Python, dicts + ndarrays, BSD, low maintenance activity), R packages (metabolomics-focused), TopSpin/nmrPipe (commercial/legacy), NMRFx/COLMARvista (recent).
- **Software design**: Three architectural decisions — MulticomplexNumbers.jl for quadrature algebra, DimensionalData.jl for named dimensions, annotation system for semantic metadata. Side-by-side code comparisons with nmrglue emphasising semantic transparency.
- **Research impact**: TITAN, NMRAnalysis.jl, NMRScreen.jl, Pluto workbooks, pH indicator database.
- **AI usage disclosure**: per JOSS 2026 requirements.

### Post-submission (JBNMR paper, later)

- NUS/IST reconstruction with multicomplex algebra.
- Annotation-driven automated processing.
- Superexperiment splitting.
- Full NMRTools + NMRAnalysis ecosystem presentation.

---

## Implementation Order

1. Repository hygiene (Phase 1) — immediate.
2. Error handling and logging audit (Phase 2) — immediate.
3. Test coverage for existing code (Phase 3) — start immediately, continue throughout.
4. Plots extension fixes (Phase 4) — early, unblocks multicomplex work.
5. `ifft!` in MulticomplexNumbers.jl — prerequisite for Phase 6.
6. Raw data loading (Phase 5) — critical dependency.
7. Processing pipeline (Phase 6) — depends on 5 and MulticomplexNumbers `ifft!`.
8. Write support (Phase 7) — overlaps with Phase 6.
9. NUS reconstruction (Phase 8) — after Phase 6.
10. Annotation extensions (Phase 9) — schema design in parallel.
11. JOSS paper (Phase 10) — after Phases 1–7.

---

## Naming Audit

The existing codebase follows Julia Base convention (lowercase, no underscores for short names). New public functions should follow suit:

| Function | Notes |
|---|---|
| `loadfid` | New. Matches `loadnmr`, `loadjdx`. |
| `savenmr` | New. Unified save, format from extension. |
| `apodize` | New public API. Wraps existing `apod`. |
| `zerofill` | New. |
| `ft` | New. Calls `fft!`/`ifft!` from MulticomplexNumbers. |
| `phase` | New. |
| `extract` | New. |
| `ist` | New. |
| `realest` | Already in MulticomplexNumbers. |
| `real`, `abs`, `imag`, `complex` | Extend `Base`. |

Internal helpers use underscore prefix: `_buildfreqdim`, `_buildtimedim`, `_replacedimtuple`, `_readfid`.