# Processing Pipeline Extensions — Proposals

Design proposals for three time-/frequency-domain processing steps to add to the composable
pipeline: **linear prediction**, **solvent suppression**, and **baseline correction**. Each follows
the conventions already established (`apodize`/`zerofill`/`ft`/`phase`): accept `NMRData` + optional
dimension identifier, return a new `NMRData`, provide curried forms for `|>`, map over
`Dict{Symbol,<:NMRData}`, keep the element type `Multicomplex` throughout, and read `:mcindex` to pick
the imaginary unit.

Recommended full pipeline ordering:

```julia
loadfid("1/ser") |>
    solventfilter |>              # time domain, on raw FID
    linearpredict(T2Dim) |>       # time domain, extend indirect dimension
    apodize |> zerofill |> ft |>
    phase(F1Dim, …) |> phase(F2Dim, …) |>
    baselinecorrect |>            # frequency domain
    real
```

## Shared infrastructure

Both LP and solvent suppression operate on the **complex** representation of one dimension's
quadrature. We already have the reinterpretation machinery used by `_mcfft!`
(`MulticomplexNumbers.unsafe_ascomplex!`/`unsafe_fromcomplex!`). Propose a small internal helper:

```julia
# run `f!(complexview)` on the Complex reinterpretation of `A` along the im-unit `unit`, in place
_withcomplex!(f!, A::AbstractArray{<:Multicomplex}, unit, dims) -> A
```

so each new operation only writes the 1-D complex kernel and lets the helper handle the multicomplex
bookkeeping (exactly as `_mcfft!` does). `frequencyshift` (below) is also shared by solvent
suppression.

A `frequencyshift(spec, dim, Δ)` primitive — multiply the FID by `exp(-imₖ·2π·Δ·t)` to move a target
frequency to zero (Δ in Hz, or resolve from a ppm shift) — is needed by solvent suppression and is
independently useful (e.g. re-referencing in the time domain). Propose exposing it publicly.

---

## 1. Linear prediction — `linearpredict`

### What others do
- **nmrPipe `LP`**: forward/backward, mirror-image; `order` (# coefficients), `pred` (# points),
  root fixing (reflect unstable poles). Common for extending truncated indirect dimensions and for
  backward-predicting corrupted first points.
- **TopSpin** (`ME_mod`/LPfr/LPbc): `NCOEF` coefficients, forward for extension, backward+complex to
  repair initial points.
- **nmrglue `lp`**: `order`, `pred`, `mode` f/b, `append` before/after, `method` svd/qr/cholesky/tls,
  `fix_roots`, forward–backward averaging, mirror-image (`lp2d`).

### Primary uses
1. **Forward extension** of a truncated dimension (usually indirect) → resolution/lineshape
   improvement without the truncation wiggles apodization would otherwise hide.
2. **Backward prediction** to replace the first few corrupted points (ringdown, baseline roll).
3. **Mirror-image LP** when the first point is exactly `t=0` (`p1=0`) or half-dwell — doubles the
   effective data.

### Proposed API

```julia
linearpredict(spec, [dim];
              npoints   = 2 * td,      # target length after prediction (forward)
              ncoef     = :auto,       # number of LP coefficients (model order)
              direction = :forward,    # :forward | :backward | :forwardbackward
              replace   = 0,           # backward: number of leading points to overwrite
              method    = :svd,        # :svd | :tls  (HSVD/LP-TLS for noisy data)
              rootfix   = :decreasing) # :decreasing | :none | :increasing
```

- `ncoef=:auto` → heuristic (e.g. `clamp(td ÷ 2, 8, 64)`); over-determined SVD is robust.
- `rootfix=:decreasing` reflects growing poles onto/into the unit circle so extrapolation can't
  blow up (the standard "predict-only-decaying-signals" stabilisation). `:increasing` for backward.
- Curried `linearpredict(dim; kw...)`, no-arg (all non-direct Fourier dimensions — LP the *direct*
  dimension is rarely wanted, so **default no-arg targets indirect Fourier dims only**), dict form.
- Updates `:td` metadata (real acquired points grow) so downstream `apodize` `tmax` stays correct.

### Multicomplex handling
LP runs on the complex (imₖ) series along the target dimension via `_withcomplex!`; the other
quadratures ride along as independent vectors. Root-fixing and coefficient solving are per-1-D-vector
(or pooled across vectors for a shared LP model — propose per-vector first, pooled as an option
later).

### Open questions / recommendation
- **Default `direction`**: `:forward`. Backward is opt-in for point repair.
- **Pooled vs per-vector coefficients**: start per-vector (simplest, matches nmrglue slice-wise);
  add `pool=true` (one AR model from all vectors, à la 2-D LP) as a follow-up.
- Ship `:svd` first; `:tls`/HSVD later.

---

## 2. Solvent suppression — `solventfilter`

### What others do
- **nmrPipe `SOL`**: time-domain. Mode 1 = **low-pass convolution difference** (subtract a moving
  average of length `-fl`, default ±16, shape `-fs` boxcar/sine/sine²); mode 2 = spline; mode 3 =
  polynomial. Removes the *on-resonance* (carrier) component, which is slowly varying in the FID.
- **TopSpin**: convolution difference; `qfil` frequency-domain notch.
- Modern: **HSVD** models the solvent as a few Lorentzians and subtracts them (frequency-selective,
  less distortion of nearby peaks).

### Principle
A signal on resonance (at the carrier) is a **DC / low-frequency** component of the FID. A low-pass
filter estimates it; subtracting leaves the off-resonance signal. For a solvent **not** at the
carrier, frequency-shift so the target sits at 0 Hz, filter, then shift back.

### Proposed API

```julia
solventfilter(spec, [dim];
              shift  = :carrier,   # :carrier | δ::ppm | (Δ::Hz)  — target frequency to suppress
              width  = 40,         # suppression notch half-width in Hz (mapped to filter length)
              shape  = :sine,      # :boxcar | :sine | :sine²  (low-pass window)
              method = :lowpass)   # :lowpass | :polynomial | :hsvd
```

- `shift=:carrier` (default) — no frequency shift needed. A ppm value (e.g. `4.7`) or an explicit Hz
  offset triggers `frequencyshift` → filter → shift back, satisfying the "specify other chemical
  shifts" requirement.
- `width` in **Hz** (intuitive: how wide a notch), converted internally to a filter length in points
  via the dimension's `swhz`. Larger `width` → shorter filter → **wider** notch; documenting the
  direction clearly matters. (Alternatively expose `filterlength` in points for nmrPipe parity — I'd
  keep both, `width` primary.)
- `shape` selects the low-pass window (boxcar = plain moving average; sine/sine² taper the notch
  edges, reducing ringing near the solvent — matches nmrPipe `-fs`).
- `method=:lowpass` first (convolution difference — simple, robust). `:hsvd` later for
  frequency-selective removal with minimal distortion of nearby peaks; `:polynomial` for parity.
- Applies to the **direct dimension by default** (solvent is in the acquisition dimension); explicit
  dim argument allowed. Time domain — before `apodize`.

### Multicomplex handling
Same `_withcomplex!` reinterpretation as LP. The frequency shift multiplies by a multicomplex phase
ramp `exp(-imₖ·2π·Δ·t)` (reuses the `phase`/`_mcfft!` unit machinery).

### Recommendation
- Default `method=:lowpass`, `shape=:sine`, `shift=:carrier`.
- Provide `frequencyshift` as a public primitive so users can also suppress multiple solvents by
  chaining, and re-reference FIDs.
- Consider `regions`/multiple targets later; start with a single target.

---

## 3. Baseline correction — `baselinecorrect`

### What others do
- **nmrPipe `POLY`/`BASE`**: polynomial (default 4th order) over automatically-detected or
  node-specified baseline regions; complex data fit real+imag with the *same* nodes. `MED` (median),
  `CBF` (constant from a region).
- **TopSpin**: `abs` auto polynomial (`ABSG` degree, `ABSF` region), spline, manual `bcm`.
- **Modern** (NMRium, metabolomics): **asymmetric least squares (ALS / airPLS)** — penalised
  iterative smoother, fully automatic, no node selection; parameters λ (smoothness), p (asymmetry),
  iterations. Rolling-ball / Whittaker also common.

### Proposed API

```julia
baselinecorrect(spec, [dim];
                method     = :als,     # :als | :polynomial | :spline | :median | :offset
                regions    = :auto,    # :auto (detect signal-free) | vector of ppm intervals
                degree     = 4,        # polynomial/spline order
                λ          = 1e6,      # ALS smoothness
                p          = 0.01,     # ALS asymmetry
                iterations = 20)
```

- **`method=:als` as the recommended default**: automatic, node-free, robust for crowded spectra —
  the main advantage over TopSpin/nmrPipe's node-based polynomial. `:polynomial` for nmrPipe/TopSpin
  parity and simple 1-D cases; `:offset` = subtract mean of a signal-free region (CBF-style).
- **`regions=:auto`**: detect baseline (signal-free) points with a noise-threshold + smoothed
  derivative classifier, reusing `estimatenoise`. `regions=[10.5 .. 11.0, -1.0 .. -0.5]` to specify
  explicit baseline windows / fit only there.
- Frequency domain — after `ft`/`phase`. Operates along the direct dimension by default.

### Multicomplex handling
Baseline is naturally a real-spectrum operation. Two options:
1. **Recommended**: require/assume the *real* part — apply after `real`, on `Float64` data, matching
   how spectroscopists baseline-correct. Simple and unambiguous.
2. **nmrPipe-style**: on `Multicomplex` data, fit the baseline to `realest` and subtract the same
   real polynomial from every component (keeps quadrature for later re-use).

Propose supporting `NMRData{<:Real}` (option 1) first; add the multicomplex variant if a use case
needs it. This is the one step where breaking the "element type never changes" rule is natural —
document that `baselinecorrect` is typically the last step before/after `real`.

### Recommendation
- Default `:als`; ship `:als` + `:polynomial` + `:offset` initially.
- Auto region detection shared with a future `:auto` phase/peak-picking noise model.

---

## Naming summary

| Function | Domain | Notes |
|---|---|---|
| `linearpredict` | time | forward/backward/mirror LP; `lp` alias optional |
| `solventfilter` | time | `shift`/`width`/`shape`/`method`; default carrier |
| `frequencyshift` | time | shared primitive, also re-referencing |
| `baselinecorrect` | frequency | ALS default; polynomial/offset/spline |

All lowercase, no underscores, consistent with `apodize`/`zerofill`/`ft`/`phase`. Internal helpers
(`_withcomplex!`, region detection, LP solvers, ALS solver) keep the underscore prefix.

## References
- nmrglue LP: <https://nmrglue.readthedocs.io/en/latest/reference/proc_lp.html>
- nmrPipe SOL: <https://spin.niddk.nih.gov/NMRPipe/doc/nmrFn/SOL.html>
- nmrPipe POLY/BASE: <https://nmrwiki.org/wiki/index.php?title=NMRPipe%3APOLY_and_BASE_baseline_correction>
- NMRium baseline (ALS/polynomial): <https://docs.nmrium.org/20_1d-spectra/baseline/>
