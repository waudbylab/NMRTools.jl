# Utilities

This page covers tools that work alongside the main data types: combining raw Bruker experiments, and importing assigned chemical shifts from the BMRB.

## Combining NMR experiments

`sumexpts(INPUTFILES...; out=OUTPUT)` is a utility function provided by NMRTools that allows you to combine multiple Bruker NMR experiments with optional weighting factors.
The function accepts the following parameters:

* `INPUTFILES...`: One or more experiment numbers or paths to be combined
* `out`: The experiment number or path where the combined data will be stored (keyword argument)

The function works with both 1D experiments (using `fid` files) and multidimensional experiments (using `ser` files).

### Basic usage

The most common way to use `sumexpts` is with experiment numbers, which is the standard Bruker convention.

```julia
julia> using NMRTools

julia> cd("myexperiment")

julia> sumexpts(11, 12, out=999)
=== NMR Experiment Summation ===
Output: 999
Inputs: 11, 12
Weights: 1.0, 1.0
Created output experiment from 11
Data format: Float64 (little-endian)
Detected nD experiment, using ser files
Cleaning processed data files in pdata/1...
  Removed: 3rrr
Removed pdata/101 directory
Removed pdata/999 directory
Loading 11/ser
Loading 12/ser
Applied weight 1.0 to experiment 11
Applied weight 1.0 to experiment 12
Writing output to 999/ser
Scans from 11: 8 × 1.0 = 8.0
Scans from 12: 16 × 1.0 = 16.0
Updated number of scans to 24
Updated title file with experiment summary
=== Operation completed successfully ===
```

This creates a new experiment directory by copying experiment 11 as a template, then adding the raw data from all input experiments together. The result is stored in experiment 999.

You can also use full paths if your experiments are in different directories:

```julia
# Combine experiments from a specific path
julia> sumexpts("myexperiment/11", "myexperiment/12", out="myexperiment/999")
# Output similar to above example
```

### Adding multiple experiments

To combine more than two experiments, provide all experiment numbers:

```julia
# Add three experiments together
julia> cd("myexperiment")

julia> sumexpts(11, 12, 13, out=999)
# Output shows all three experiments being combined
```

By default, all experiments are weighted equally. The result will be the sum of all input data.

### Creating difference spectra

To create a difference spectrum, use the `weights` parameter:

```julia
# Create a difference spectrum (experiment 12 subtracted from experiment 11)
julia> cd("myexperiment")

julia> sumexpts(11, 12, out=999, weights=[1.0, -1.0])
=== NMR Experiment Summation ===
Output: 999
Inputs: 11, 12
Weights: 1.0, -1.0
Created output experiment from 11
Data format: Float64 (little-endian)
Detected nD experiment, using ser files
Cleaning processed data files in pdata/1...
  Removed: 2rr
Loading 11/ser
Loading 12/ser
Applied weight 1.0 to experiment 11
Applied weight -1.0 to experiment 12
Writing output to 999/ser
Scans from 11: 16 × 1.0 = 16.0
Scans from 12: 16 × 1.0 = 16.0
Updated number of scans to 32
Updated title file with experiment summary
=== Operation completed successfully ===
```

This is particularly useful for applications like:
- Reaction monitoring
- Comparing before/after spectra
- Highlighting changes between samples
- Removing background signals

### Weighted combinations

You can apply different weights to each experiment to create weighted averages or to emphasize certain spectra:

```julia
# Weighted combination of three experiments
julia> cd("myexperiment")

julia> sumexpts(11, 12, 13, out=999, weights=[0.5, 1.0, 2.0])
=== NMR Experiment Summation ===
Output: 999
Inputs: 11, 12, 13
Weights: 0.5, 1.0, 2.0
Created output experiment from 11
Data format: Float64 (little-endian)
Detected nD experiment, using ser files
Cleaning processed data files in pdata/1...
Loading 11/ser
Loading 12/ser
Loading 13/ser
Applied weight 0.5 to experiment 11
Applied weight 1.0 to experiment 12
Applied weight 2.0 to experiment 13
Writing output to 999/ser
Scans from 11: 8 × 0.5 = 4.0
Scans from 12: 16 × 1.0 = 16.0
Scans from 13: 32 × 2.0 = 64.0
Updated number of scans to 84
Updated title file with experiment summary
=== Operation completed successfully ===
```

The weights vector must have the same length as the number of input experiments.

### Experimental details

The function handles several important details automatically:

- The number of scans (`NS` parameter) in the output is updated to the weighted sum of input scan counts
- The title file is updated to indicate which experiments were combined and their weights
- Processed data files in `pdata/1` are removed to prevent inconsistencies
- If experiments have different data sizes, they are truncated to match the smallest dataset

### Command-line usage

The function can also be called from the command line:

```bash
cd myexperiment
julia -e 'using NMRTools; sumexpts(11, 12, out=999)'
```

Or for creating a difference spectrum:

```bash
cd myexperiment
julia -e 'using NMRTools; sumexpts(11, 12, out=999, weights=[1.0, -1.0])'
```

This allows for easy integration with processing workflows and shell scripts.

## Assigned chemical shifts from the BMRB

`loadbmrb` downloads the assigned chemical shifts of a [BMRB](https://bmrb.io) entry and
assembles them into the peaks expected in a given experiment, as a peak list of the kind described
in [Peaks and peak lists](tutorials/peaks.md):

```julia
hsqc = loadbmrb(15000, H1, N15)          # ¹H,¹⁵N HSQC, side chains included
nmrplot(hsqc; xlabel="¹H (ppm)", ylabel="¹⁵N (ppm)")

loadbmrb(15000, H1, C13)                 # ¹H,¹³C HSQC
```

Called with just an entry ID (or the name of a local NMR-STAR file), `loadbmrb` returns the
shifts themselves. Download them once and build several lists from them with `peaklist`,
which takes the same dimension arguments:

```julia
shifts = loadbmrb(15000)
hsqc = peaklist(shifts, "H", "N")        # backbone amide ¹H,¹⁵N HSQC
```

Dimensions are given in spectrum order, direct dimension first. Each argument names an atom,
optionally with an offset to a neighbouring residue, so any triple-resonance correlation can
be written down directly:

```julia
peaklist(shifts, "CO-1", "N")            # CON: C′ of residue i-1 with N of residue i
peaklist(shifts, "H", "N", "CO-1")       # HNCO
peaklist(shifts, "H", "N", "CA")         # HNCA, i peaks
peaklist(shifts, "H", "N", "CA-1")       # HNCA, i-1 peaks
```

Passing nuclei instead of atom names gives every directly bonded X–H pair, which is what an
HSQC actually shows. This includes side-chain NH and NH₂ groups for ¹⁵N, and every CH, CH₂
and CH₃ group for ¹³C. Methyl and methylene protons with identical shifts give a single peak:

```julia
peaklist(shifts, H1, N15)                # ¹H,¹⁵N HSQC including side chains
peaklist(shifts, H1, C13)                # ¹H,¹³C HSQC
```

A shift list is a `Vector{ChemicalShift}`, so residues can be selected with `filter` before
building peaks, e.g. for the methyl region of an ILV-labelled sample:

```julia
ilv = filter(s -> s.residuename in ("ILE", "LEU", "VAL"), shifts)
peaklist(ilv, H1, C13)
```

Peaks are labelled by residue and atom (`L23H-N`, `D5C-F6N`, `V9HG2-CG2`), so a predicted
list can be overlaid on a spectrum with `nmrplot!`, or saved with `savepeaks` as the starting
point for an analysis in NMRAnalysis.jl. BMRB shifts were measured under the deposited sample
conditions and referencing, so expect systematic offsets; `register!` finds and removes the
overall offset (see [Registering a peak list onto a spectrum](@ref)):

```julia
pks = loadbmrb(15000, H1, N15)
register!(pks, spec)
```
