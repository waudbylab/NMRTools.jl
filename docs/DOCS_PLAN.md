# Documentation restructure plan

Status: implemented on branch `docs/restructure`. This file lives outside `docs/src` so it is not
part of the built site. Differences from the plan as written below:

- The BMRB 17769 file is in the repository (`data/bmr17769_3.str`), not an artifact.
- The API has an extra Plotting page (`nmrplot`, `nmrplot!`).
- The BMRB tutorial assigns the ¹H,¹⁵N projection of `3D_HNCA` (298 K, matching the BMRB entry);
  `2D_HN` (277 K) is used to show the limits of a uniform offset.
- The solvent filter and 2D autophasing problems found while writing the examples were fixed in the
  package (see `NEWS.md`), so the examples use `autoprocess` with its defaults. The 2D examples
  still show TopSpin's phases where an exact result matters.

Baseline: the docs on `main` at v0.3.0 (`66e6fbd`). This feature branch's `docs/` is behind `main`
(main has deleted `tutorials/metadata.md` and revised `peaks.md`, `processing.md` and `utilities.md`),
so the plan is written against `main`.

## Decisions so far

- Add `exampledata(name; fid=true)` so raw-data examples run during the build.
- `sumexpts` and BMRB peak lists get a top-level **Tools** section with executed examples, linked
  from Home and Getting started.
- Bundle the NMR-STAR file for ubiquitin, BMRB entry 17769, and add a BMRB tutorial on it. No
  synthetic spectra.
- Makie for every plot, except on the Plots.jl page.
- No Worked examples section. Relaxation stays as a tutorial; diffusion and titration are dropped
  (NMRAnalysis.jl covers proper analysis).
- Getting started stays a top-level item.
- Tutorials on working with data and on processing plot the data at each stage (FIDs, 1Ds, 2Ds).
- Peak detection, registration and peak list I/O get how-to guides; registration is taught in the
  BMRB tutorial only.
- "Reference chemical shifts", not "Reference spectra".
- Annotations go under an **Experimental** heading.
- No redirects for old URLs.

## How to decide what kind of page something is

| Kind | Reader | Shape | Title | Test |
|---|---|---|---|---|
| Tutorial | new to NMRTools | one dataset, one path, default settings, read in order, says why as it goes | the activity: "Processing a 1D spectrum" | the default route |
| How-to | knows the basics, has a task | one task with its options and variants, any order, little explanation, links to concepts | a verb: "Load spectra", "Plot 2D spectra" | a knob on the default route |
| Reference | looking something up | complete tables and lists, not read through | the thing: "Window functions" | the list of every knob |
| Concept | wants to understand | explanation, no tasks | the idea: "Dimensions" | why the knob exists |

Applied to processing: `autoprocess` on a 1D is the tutorial; manual step-by-step processing and a 2D
with linear prediction are how-tos; the catalogue of window functions is reference; the digital
filter and group delay are a concept.

The Tools section is the one deliberate exception: it is there for visibility. Its pages are
how-to-shaped (all the options of one tool).

## Example data

| Dataset | Contents | Raw data | Used by |
|---|---|---|---|
| `1D_1H` | ubiquitin, 10% D₂O | `fid` | Getting started, tutorial 2 |
| `1D_19F` | ¹⁹F 1D | `fid` | Getting started, tutorial 1 |
| `1D_19F_titration` | 11 × ¹⁹F 1D | 11 × `fid` | Combining experiments, Animate plots |
| `2D_HN` | ubiquitin ¹H,¹⁵N | `ser` | Getting started, tutorials 1, 3, 4, several how-tos |
| `2D_HN_titration` | series of ¹H,¹⁵N | — | Plot 2D spectra (overlays), Animate plots |
| `3D_HNCA` | ubiquitin BEST-HNCA, 298 K | `ser` | tutorial 4, Load spectra |
| `3D_HNCA_pipe`, `3D_HNCA_ucsf` | the same, nmrPipe and UCSF | — | Load spectra |
| `pseudo2D_T2` | relaxation series | `ser` | tutorials 1 and 5, Process step by step |
| `pseudo2D_XSTE`, `pseudo3D_*` | diffusion, R2, kinetics | some | Plot 2D spectra (pseudo-2D) |
| BMRB 17769 (new) | ubiquitin shifts | — | tutorial 4, Peak lists from the BMRB |

BMRB 17769 (Cornilescu, Marquardt, Ottiger, Bax) is wild-type human ubiquitin, ¹³C,¹⁵N-labelled,
pH 6.5, 298 K, 50 mM phosphate, 7% D₂O, with 570 ¹H, 348 ¹³C and 83 ¹⁵N shifts. It is the deposition
behind the 1D3Z structure and matches the `3D_HNCA` temperature. (5387, the other complete entry,
is in reverse micelles.)

## Proposed navigation

```
Home                                  index.md
Getting started                       getting-started.md

Tutorials
  1 Working with spectra              tutorials/spectra.md
  2 Processing a 1D spectrum          tutorials/processing.md
  3 Peaks and peak lists              tutorials/peaks.md
  4 Assigning ubiquitin from the BMRB tutorials/bmrb.md
  5 Relaxation analysis               tutorials/relaxation.md

Tools
  Combining experiments               tools/sumexpts.md
  Peak lists from the BMRB            tools/bmrb.md

How-to guides
  Data
    Load spectra                      howto/load.md
    Process step by step              howto/process.md
    Process a 2D with linear prediction  howto/linearprediction.md
    Reference chemical shifts         howto/referencing.md
    Save spectra and FIDs             howto/save.md
  Plotting
    Plot 1D spectra                   howto/plot1d.md
    Plot 2D spectra                   howto/plot2d.md
    Animate plots                     howto/animate.md
    Plot with Plots.jl                howto/plots.md
  Peaks
    Detect peaks                      howto/detectpeaks.md
    Register a peak list onto a spectrum  howto/register.md
    Save and load peak lists          howto/peaklists.md
  Metadata
    Find acquisition parameters       howto/acqus.md
    Find sample information           howto/samples.md

Concepts
  The NMRData model                   concepts/nmrdata.md
  Dimensions                          concepts/dimensions.md
  Multicomplex numbers and quadrature concepts/multicomplex.md
  Metadata                            concepts/metadata.md
  Time and frequency units            concepts/units.md
  Bruker data and TopSpin conventions concepts/bruker.md

Reference
  Metadata keys                       reference/metadatakeys.md
  Window functions                    reference/windows.md
  Power levels                        reference/power.md
  Frequency lists                     reference/fqlists.md
  Sample metadata schema              reference/samples.md
  Nuclei, coherences and Xi ratios    reference/nuclei.md
  Makie plot options                  reference/makieoptions.md

Experimental
  Pulse programme annotations         experimental/annotations.md

API
  Loading and saving                  api/io.md
  Processing                          api/processing.md
  Peaks and assignments               api/peaks.md
  Dimensions                          api/dimensions.md
  Metadata                            api/metadata.md
  Window functions and lineshapes     api/windows.md
  Nuclei and coherences               api/nuclei.md
  Internals                           api/internals.md

Release notes                         releasenotes.md
```

The how-to subgroups follow the order work is done in: get data in, turn it into a spectrum,
reference it, save it; then plot it; then pick and assign peaks; then look things up. Documenter
supports nested page groups.

## Page contents

### Getting started

Ten minutes, Makie throughout: install; load and plot `1D_19F`; zoom with `..`; load and plot
`2D_HN`; process `1D_1H` from its FID in one line with `autoprocess`; links to the tutorials and both
Tools pages. Don't print the full `NMRData` (it dumps 19 metadata entries and the CI runner's
artifact path); show selected fields.

### Tutorials

Tutorials 1 and 2 plot the data at each stage with `nmrplot`, so plotting is learnt along the way
and the plotting how-tos only add options.

1. **Working with spectra** (`1D_19F`, `2D_HN`, `pseudo2D_T2`): what an `NMRData` is, array and
   selector indexing, axes and `data`, real and imaginary parts, a first look at metadata, saving in
   one line. Plots: the 1D spectrum and a selected region; the 2D spectrum and a sub-region; a 1D
   slice taken from the 2D; real and imaginary parts of the 1D; a pseudo-2D series. Links to the
   concepts pages.
2. **Processing a 1D spectrum** (`1D_1H` FID): `autoprocess`, then what it did. Plots: the raw FID
   (including the digital-filter build-up before t = 0); the FID after the solvent filter; the
   window function over the FID; the zero-filled FID; the spectrum straight after `ft`; the phased
   spectrum. Ends by linking to the processing how-tos.
3. **Peaks and peak lists** (`2D_HN`): `detectpeaks`, the peak list as a table, filtering, plotting
   peaks on the spectrum, `findpeak`, saving and loading. Registration is left to tutorial 4.
4. **Assigning ubiquitin from the BMRB** (`2D_HN`, `3D_HNCA`, BMRB 17769):
   1. load the shift list and show it as a table
   2. predict the ¹H,¹⁵N HSQC and overlay it on `2D_HN`; note the systematic offset
   3. `register!` the prediction onto the spectrum; print the fitted offset; replot
   4. use the registered labels as assignments; zoom to a few residues
   5. side chains: backbone-only (`"H", "N"`) against all NH groups (`H1, N15`), showing the Asn/Gln
      NH₂ and Arg/Trp side-chain peaks
   6. check assignments in the HNCA: `"CA"` and `"CA-1"` predictions on ¹³C strips from `3D_HNCA`
   7. save the assigned peak list with `savepeaks` for NMRAnalysis.jl
5. **Relaxation analysis** (`pseudo2D_T2`): integrate across the series, fit with LsqFit, plot. Says
   at the top that NMRAnalysis.jl does this properly; this page shows how the pieces fit together.

### Tools

- **Combining experiments** (copy of `1D_19F_titration` in a temporary directory): sum two
  experiments into a new experiment number with the summary output shown once; process and overlay
  the result; difference spectrum (weights `1, -1`); weighted combination of three; command-line
  usage (not executed).
- **Peak lists from the BMRB** (BMRB 17769): every way of writing a correlation, each executed and
  plotted: `loadbmrb(file, H1, N15)`, loading shifts once and building several lists with
  `peaklist`, atom offsets (CON `"CO-1", "N"`, HNCO, HNCA i and i−1), nucleus arguments for every
  bonded X–H pair, ¹H,¹³C HSQC, an ILV-methyl subset with `filter`, peak labels, `register!`,
  `savepeaks`. The download form `loadbmrb(17769, ...)` as a non-executed block at the top.

The tutorial tells the ubiquitin story once; the Tools page is the complete set of options. They
share the data file and the tutorial links to the Tools page for anything it doesn't cover.

### How-to guides

**Data**

- **Load spectra**: Bruker processed data (experiment, pdata), raw data with `loadfid` and what
  happens on loading, nmrPipe (including plane series with `%03d` templates), UCSF, the
  `experimentfolder` keyword. (`2D_HN`, `3D_HNCA`, `3D_HNCA_pipe`, `3D_HNCA_ucsf`)
- **Process step by step**: `apodize`, `zerofill`, `ft`, `phase` and `solventfilter` by hand;
  choosing a window function; manual and automatic phasing; `extract`; per-dimension arguments;
  a pseudo-2D series, where only the Fourier dimension is processed and the non-Fourier axis is set
  with `setrelaxtimes`, `setgradientlist` or `setkinetictimes`. (`1D_1H`, `2D_HN`, `pseudo2D_T2`
  FIDs)
- **Process a 2D with linear prediction**: `autoprocess` with per-dimension `Dict` settings,
  `linearpredict` in the indirect dimension, comparison with and without. (`2D_HN` ser)
- **Reference chemical shifts**: `old => new`, indirect referencing with Xi ratios, DSS vs TMS,
  water referencing from the sample temperature. (`2D_HN`)
- **Save spectra and FIDs**: `savenmr` and `savefid`, `procno`, what is copied from the source
  experiment, opening the result in TopSpin, the sign convention in one paragraph.

**Plotting** (Makie except the last)

- **Plot 1D spectra**: limits, overlays, normalisation, stacked plots, labels, saving, multi-panel
  layouts.
- **Plot 2D spectra**: contour levels and colours, overlays, normalisation, pseudo-2D data, GLMakie
  keyboard navigation, saving, layouts.
- **Animate plots**: `record` over `2D_HN_titration` and `1D_19F_titration`.
- **Plot with Plots.jl**: the Plots.jl recipes, 1D and 2D, with their options table.

**Peaks**

- **Detect peaks**: `detectpeaks` thresholds (in units of noise), positive and negative peaks,
  restricting to a region, 1D, 2D and 3D data, `findpeak` for a single peak near a position,
  `findnearest`, filtering a peak list by intensity, linewidth or label. (`1D_19F`, `2D_HN`,
  `3D_HNCA`)
- **Register a peak list onto a spectrum**: `register!` and `register` with a peak list from the
  BMRB, from a saved file or from another spectrum; the `scale` argument; checking the fitted
  offset; registering one list across a series of spectra. (`2D_HN`, `2D_HN_titration`, BMRB 17769)
- **Save and load peak lists**: `savepeaks` and `loadpeaks`, NMRTools CSV, NMRAnalysis.jl and Sparky
  formats, peak lists as Tables.jl tables (e.g. into a `DataFrame`).

**Metadata**

- **Find acquisition parameters**: `acqus`, arrayed parameters, `channels`, power levels, frequency
  lists, vdlists. (`2D_HN`)
- **Find sample information**: `sample` and its accessors, `scanexperiments`, `scansamples`,
  `findsample`, `findexperiments`.

### Concepts

- **The NMRData model**: `NMRData` as an `AbstractDimArray`; what operations preserve; `rebuild`;
  link to the DimensionalData.jl docs.
- **Dimensions**: the dimension type tree (diagram), F/T/G/X dimensions, how `ft` turns a `T`
  dimension into an `F` dimension, `replacedimension`, `finddim`/`resolvedim`.
- **Multicomplex numbers and quadrature**: one imaginary unit per dimension; `im1` not `im`;
  component ordering; `:mcindex`; the element type stays the same through processing; `real`,
  `realest`, `complex`. Uses `multicomplex.png`.
- **Metadata**: spectrum vs axis metadata (explained once here), lookup syntax, the special metadata
  types with links to their reference pages.
- **Time and frequency units**: all times in seconds and frequencies in Hz (pulse lengths, delays,
  vdlists, `bf`/`sf`), where Bruker uses µs, ms and MHz; `ppm` and `hz` conversions.
- **Bruker data and TopSpin conventions**: digital filter and group delay, baseopt, `DSPFVS` edge
  cases, `2^NC` and first-point scaling, quadrature modes (States, States-TPPI, echo/anti-echo), sign
  convention and why `PHC0`/`PHC1` carry over, `allcomponents=true`. Some of this is currently only
  in `CLAUDE.md`.

### Reference

Metadata keys (tables of standard spectrum and axis keys), window functions (each type with formula
and plot), power levels, frequency lists, the sample metadata schema, nuclei and coherences with the
Xi ratio tables, and the Makie `nmrplot` keyword tables.

### Experimental

**Pulse programme annotations**: an introduction and one worked annotated pulse programme, then the
schema (v0.0.2).

## Where current pages go

| Current page (main) | New location |
|---|---|
| `index.md` | rewritten: figure, runnable Makie example, feature list linking to Tools |
| `quickstart.md` | `getting-started.md` (rewritten) |
| `tutorials/data.md` | `tutorials/spectra.md`; multicomplex material to `concepts/multicomplex.md`; loading to `howto/load.md`; referencing to `howto/referencing.md`; saving to `howto/save.md` |
| `tutorials/processing.md` | `autoprocess` to `tutorials/processing.md`; individual steps and pseudo-2D series to `howto/process.md`; per-dimension settings to `howto/linearprediction.md`; digital filter and scaling to `concepts/bruker.md` |
| `tutorials/makie.md` | `howto/plot1d.md`, `howto/plot2d.md`; options tables to `reference/makieoptions.md`; animations to `howto/animate.md` |
| `tutorials/plotting.md` | `howto/plots.md` |
| `tutorials/peaks.md` | `tutorials/peaks.md`; detection options to `howto/detectpeaks.md`; saving and loading formats to `howto/peaklists.md`; registration to `tutorials/bmrb.md` and `howto/register.md` |
| `tutorials/relaxation.md` | `tutorials/relaxation.md` (Makie plots) |
| `tutorials/diffusion.md` | deleted |
| `tutorials/animation.md` | `howto/animate.md` |
| `metadata/index.md` | concept prose to `concepts/metadata.md`; acqus, channels, aux files to `howto/acqus.md` |
| `metadata/timeunits.md` | `concepts/units.md` ("Time and frequency units") |
| `metadata/power.md` | `reference/power.md` |
| `metadata/fqlists.md` | `reference/fqlists.md` (plus "Frequency list API" from `reference/metadata.md`) |
| `metadata/windows.md` | `reference/windows.md` |
| `metadata/samples.md` | schema to `reference/samples.md`; accessors and scanning to `howto/samples.md` |
| `utilities.md` | `tools/sumexpts.md`, `tools/bmrb.md` (both rewritten with executed examples) |
| `reference/referencing.md` | `howto/referencing.md`; Xi ratio tables to `reference/nuclei.md` |
| `reference/metadata.md` | tables to `reference/metadatakeys.md`; prose to `concepts/metadata.md`; `@docs canonical=false` blocks dropped |
| `reference/annotations.md` | `experimental/annotations.md` |
| `reference/coherences.md` | `reference/nuclei.md` |
| `api.md` | split into `api/*.md` |
| `indexes.md` | deleted |

## API

One page per topic using `@docs` lists, with an `@index` block at the top of each page. Remove the
`size_threshold` override in `make.jl`. Use `checkdocs=:exports` so that an exported name missing
from every page fails the build, instead of the current `warnonly=[:missing_docs]`.

- **Loading and saving:** `loadnmr`, `loadfid`, `savenmr`, `savefid`, `loadmetadata`, `loadjdx`, `exampledata`, `sumexpts`, `withreference`
- **Processing:** `apodize`, `zerofill`, `ft`, `phase`, `extract`, `frequencyshift`, `linearpredict`, `solventfilter`, `autoprocess`, `Pole`, `estimatepoles`, `estimatephase`, `autophase`, `autophase0`, `estimatenoise`, `decimate`, `scale`, `reference`, `isaqueous`, `watershift`
- **Peaks and assignments:** `Peak`, `intensity`, `linewidth`, `detectpeaks`, `findpeak`, `findnearest`, `register`, `register!`, `peaklist`, `ChemicalShift`, `loadpeaks`, `savepeaks`, `loadbmrb`, `pcls`
- **Dimensions:** `NMRData`, `AbstractNMRData`, `NMRDimension` and subtypes, `F1Dim`…`SpinlockDim`, `ppm`, `hz`, `replacedimension`, `set*` functions, `shiftdim`, `finddim`, `resolvedim`, `hasnonfrequencydimension`
- **Metadata:** `metadata`, `defaultmetadata`, `metadatahelp`, `label`, `label!`, `units`, `acqus`, `channels`, `channel`, `annotations`, `hasannotations`, `sample`, `hassample`, `samplefile`, `referencepulse`, `NMRSample`, `NMRExperiment`, `FQList`, `Power`, `db`, `watts`, `scanexperiments`, `scansamples`, `findsample`, `findexperiments`
- **Window functions and lineshapes:** `WindowFunction` and subtypes, `apod`, `lineshape`, `LineshapeComplexity`, `RealLineshape`, `ComplexLineshape`
- **Nuclei and coherences:** `Nucleus`, `nucleus`, `spin`, `gyromagneticratio`, `xi_ratio`, `Coherence`, `SQ`, `MQ`, `coherenceorder`
- **Internals:** documented non-exported functions and `SchemaMigrate`, marked as not public API.

DimensionalData re-exports (`At`, `Near`, `..`, `dims`, `rebuild`, …) are not documented here; the
Dimensions page links to the DimensionalData.jl docs.

## Doctests

Documenter runs `jldoctest` blocks in docstrings during `makedocs`. There are five today, in
`coherences.jl` (`MQ`, `coherenceorder`) and `isotopedata.jl` (`spin`, `gyromagneticratio`,
`xi_ratio`). Most of the other 46 docstring examples are plain `julia` or `julia-repl` blocks, and
several of the `julia-repl` ones show output that is never checked.

1. **Convert examples that need no data** to `jldoctest`:
   - `nucleus` (`nuclei.jl`): `nucleus("19F")` → `F19`, and the `ArgumentError` case
   - `Power`, `db`, `watts` (`power.jl`): dB ↔ W conversions, printed rounded so float formatting
     can't break them
   - `Peak` (`peaks.jl`): constructor and field access, if `show` output is stable
2. **Convert examples that read one value from example data** to `jldoctest`, with a shared setup
   that loads `exampledata("2D_HN")`:
   - `acqus(spec, :pulprog)` and an arrayed parameter (`metadata.jl`)
   - `channels(spec, C13)[:bf]` (`metadata.jl`), rounded
   - `finddim(spec, H1)`, `finddim(spec, N15)` (`dimensions.jl`); this example currently uses a
     hypothetical `hsqc.ft2`
   - `referencepulse` (`metadata.jl`), if the example data has the calibration it needs
3. **Leave as plain `julia` blocks**: processing pipelines, `loadfid`/`autoprocess`, peak detection,
   `loadbmrb`, `scanexperiments`, and anything that prints an `NMRData`. Update their paths to use
   `exampledata` where possible and link each to the page where it runs.
4. **Run doctests in the test suite** with `Documenter.doctest(NMRTools; manual=false)` in a
   SafeTestset, so a broken docstring example fails `Pkg.test()` and not only the docs build.

## Code prerequisites

1. **`exampledata(name; fid=true)`**: returns `loadfid` of the dataset's raw data (table above).
   `sumexpts` needs the experiment directory rather than loaded data; the `:experimentfolder`
   metadata of the result gives it, or a small `examplepath(name)` could be added.
2. **BMRB 17769 NMR-STAR file**, as an artifact (preferred, alongside the other example data) or a
   docs asset, with an `exampledata`-style accessor so docs and tests can find it.
3. **A compact `show` for `NMRData`** (optional): the default display prints every metadata entry,
   including absolute file paths.
4. **Plotting time-domain data** (check first). The processing tutorial plots FIDs at each step.
   `nmrplot` already handles a 1D FID (a 1D `NonFrequencyDimension` in `ext/MakieExt/plot_1d.jl`),
   but plots only the real component (`realest`). Showing real and imaginary parts together, and
   plotting one increment of a 2D `ser`, may need small additions. Phasing plots in the processing
   how-to want the same for spectra.

## Phases

1. **Structure.** Move and rename files, update `make.jl` `pages`, fix `@ref` links.
2. **Data.** `exampledata(name; fid=true)`, the BMRB 17769 artifact.
3. **Tools and BMRB tutorial.** `tools/sumexpts.md`, `tools/bmrb.md`, `tutorials/bmrb.md`; Home figure.
4. **Processing.** `tutorials/processing.md`, the three processing how-tos, `concepts/bruker.md`.
5. **Remaining how-tos and concepts.**
6. **API split and doctests.**
7. **DocumenterVitepress (optional)**, ideally with MulticomplexNumbers.jl, NMRAnalysis.jl and
   NMRScreen.jl at the same time:
   - `make.jl`: `format = DocumenterVitepress.MarkdownVitepress(repo=..., devbranch="main", devurl="dev")`
   - `docs/src/.vitepress/config.mts`: top nav (Getting started, Tutorials, Tools, How-to, Concepts,
     Reference, API), sidebar per section, base path for `waudbylab.org/NMRTools.jl`
   - `index.md`: hero frontmatter and feature grid
   - check `@ref` resolution, LaTeX in docstrings and the GIF/SVG outputs
   - CI: confirm build time and the versions dropdown

## Open questions

None outstanding. Resolved in the last review: tutorials 1 and 2 plot along the way; pseudo-2D
processing is part of "Process step by step"; registration is taught in the BMRB tutorial and has
its own how-to.
