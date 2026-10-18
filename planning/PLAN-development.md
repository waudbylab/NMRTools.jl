# NMRTools.jl Development Plan

## Design Principles

1. **Domain modelling first.** Types encode the structure of NMR data: multicomplex algebra for quadrature, named dimensions for physical axes, metadata for acquisition parameters. Processing functions dispatch on these types.

2. **Julian idioms.** Extend `Base` and standard library functions where appropriate (`real`, `complex`, `abs`, `fft!`). Public function names are lowercase without underscores, following Julia Base convention (`loadfid`, `apodize`, `zerofill`, `lineshape`). Types are CamelCase. Underscore-prefixed names are internal only. Use Julia standard exception types where applicable; reserve custom exceptions for domain-specific errors.

3. **Multicomplex numbers are the data representation.** `im1` is the conventional complex unit. Element types never change during processing — only `real(spec)` at the end of a pipeline strips imaginary components. There is no mixing with Julia's built-in `Complex`; use `complex(spec)` to convert a first-order multicomplex to `Complex` when needed.

4. **`:mcindex` metadata links dimensions to imaginary units.** Set once by `loadfid`, never modified, carried through all operations. Every processing function reads `:mcindex` from the target dimension's metadata to determine which imaginary unit to operate on.

5. **Composable pipeline.** Processing functions take and return `NMRData`. Curried forms enable `|>` piping. No-argument forms apply sensible defaults to all applicable dimensions. Explicit dimension arguments available when control is needed.

---

## Status at v0.3.0 (October 2026)

| Phase | Status |
|---|---|
| 1. Repository hygiene | Done. CONTRIBUTING, code of conduct, CI and Codecov badges, CITATION.cff current. |
| 2. Error handling and logging audit | Done. Conventions recorded in `CLAUDE.md`. One silent `catch` remains in `channels` (`NMRBase/metadata.jl`). |
| 3. Test coverage | Partly done. Gaps listed below. |
| 4. Plots extension fixes | Mostly done (colours, normalisation, multicomplex data, typed errors). Pseudo-2D items remain. |
| 5. Raw time-domain loading | Done for 1D to 3D: `loadfid` with quadrature recombination, group delay handling and first-point scaling. NUS, 4D and the legacy filter table remain. |
| 6. Processing pipeline | Done: `apodize`, `zerofill`, `ft`, `phase`, `extract`, `frequencyshift`, `linearpredict`, `solventfilter`, `autophase`, `autoprocess`. |
| 7. Write support | Bruker done (`savenmr` to pdata, `savefid` to fid/ser). nmrPipe and JLD2 remain. |
| 8. NUS reconstruction | Not started. Depends on NUS loading. |
| 9. Annotation schema extensions | Not started. |
| 10. JOSS paper | Not started. Remaining blocker from the minimum viable list: nmrPipe output. |

Also delivered in 0.3, beyond the original plan: the peak framework (`Peak{N}`, N-dimensional
`detectpeaks`, `register!`, peak list I/O) and BMRB chemical shift import (`loadbmrb`, `peaklist`).

The detailed designs of the completed phases have been removed from this plan. The code, its
docstrings and the processing tutorial are now the reference for them.

---

## Outstanding work

### Test coverage

- `sumexpts`: file I/O, weighted summation, metadata updates.
- Window functions: every `apod` method, including `GaussWindow`, `LorentzToGaussWindow` and `GeneralSineWindow`.
- Format readers: malformed acqus, missing parameters, unusual dimensionalities.
- `loadbmrb` against the live BMRB (issue #78).

Target: >80% line coverage, with particular depth on multicomplex and annotation paths.

### Plots extension

#### Pseudo-2D recipe

The default `seriestype --> :path3d` waterfall is unusual. Consider changing the default to a heatmap or stacked 1D overlay, which are more standard for NMR pseudo-2D data. Keep 3D waterfall available via explicit `seriestype=:path3d`.

The `usegradient` attribute is read but never works (deleted from attributes before the dispatch path reaches it). Either fix the dispatch or remove the dead code.

### Raw time-domain loading

- **NUS:** if `FnTYPE` indicates NUS or a `nuslist` file is present, load the schedule, place the
  acquired FIDs at their positions in a zero-filled grid, and record `:nus => true` and `:nuslist`
  in the dimension metadata.
- **4D data** in `loadfid`, `savefid` and `savenmr` (currently limited to 3D).
- **Legacy digital filter:** a `(DSPFVS, DECIM)` to group delay lookup table for data without
  `GRPDLY`, written from scratch citing Westler and Abildgaard. At present `loadfid` warns and
  removes no delay.

### Write support

#### Unified `savenmr`

`savenmr(path, spec; procno=1)` currently writes Bruker pdata only. Extend it, keeping that
behaviour for paths without a recognised extension, to a single save function that infers format from file extension:

```julia
savenmr("output.ft2", spec)      # nmrPipe format (from .ft1/.ft2/.ft3 extension)
savenmr("output.nmr", spec)      # JLD2 serialisation (from .nmr/.jld2 extension)
savenmr("output/pdata/1", spec)  # Bruker pdata (from path structure)
savenmr("output.ucsf", spec)     # UCSF/Sparky (future)
```

Optional keyword to override format detection: `savenmr("output.dat", spec; format=:nmrpipe)`.

Corresponding loading already handled by `loadnmr` with `getformat` dispatcher — extend `getformat` to recognise `.nmr`/`.jld2` for serialised data.

#### nmrPipe write

Populate header from NMRData metadata (the header specification is already fully parsed in `parsenmrpipeheader`). Write header + Float32 data. Support 1D, 2D (single file), 3D (plane series).

#### JLD2 serialisation

Workflow persistence for Julia-to-Julia use. Store schema version for forward compatibility. Not interoperable with other software.

Tests: nmrPipe round trip (write, read back, compare) and JLD2 round trip for each `NMRData` type.

### NUS reconstruction

#### IST

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

#### Tests

- Retrospectively undersampled 2D, compare with fully sampled.
- 100% sampling reproduces standard FT.
- Schedule loading and mask generation.

### Annotation schema extensions

#### Quadrature annotations

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

#### Processing hints (forward-looking)

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

Design the schema now. `autoprocess` exists, so it can take these hints as its defaults once the
schema is settled.

#### Phase cycle annotations (forward-looking)

```yaml
;@ phasecycle:
;@   steps: 4
;@   receiver: [0, 180, 0, 180]
;@   experiments:
;@     zqc: {weights: [1, -1, 1, -1]}
;@     dqc: {weights: [1, 1, -1, -1]}
```

Design the schema now. Implement `loadfid(; split=true)` later.

### JOSS paper

Minimum viable for submission: everything in the status table marked done, plus nmrPipe output and
test coverage above 80%.

#### Paper structure

- **Summary**: Julia library encoding the algebraic and dimensional structure of NMR data in the type system.
- **Statement of need**: NMR data has three layers of structure (multicomplex quadrature, named physical dimensions, acquisition metadata) that existing tools represent only by convention. NMRTools.jl makes them explicit.
- **State of the field**: nmrglue (Python, dicts + ndarrays, BSD, low maintenance activity), R packages (metabolomics-focused), TopSpin/nmrPipe (commercial/legacy), NMRFx/COLMARvista (recent).
- **Software design**: Three architectural decisions — MulticomplexNumbers.jl for quadrature algebra, DimensionalData.jl for named dimensions, annotation system for semantic metadata. Side-by-side code comparisons with nmrglue emphasising semantic transparency.
- **Research impact**: TITAN, NMRAnalysis.jl, NMRScreen.jl, Pluto workbooks, pH indicator database.
- **AI usage disclosure**: per JOSS 2026 requirements.

#### Post-submission (JBNMR paper, later)

- NUS/IST reconstruction with multicomplex algebra.
- Annotation-driven automated processing.
- Superexperiment splitting.
- Full NMRTools + NMRAnalysis ecosystem presentation.

## Implementation order

1. Test coverage gaps, continuing throughout.
2. nmrPipe write, which unblocks the JOSS submission.
3. NUS loading, then IST reconstruction.
4. JLD2 serialisation, 4D support and the legacy filter table, as needed.
5. Annotation extensions, with schema design in parallel.
6. JOSS paper.
