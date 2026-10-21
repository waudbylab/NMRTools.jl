# Combining experiments

[`sumexpts`](@ref NMRTools.NMRIO.sumexpts) adds raw Bruker experiments together, with optional
weights, and writes the result as a new experiment. Use it to add repeated acquisitions, to make
difference spectra, or to combine experiments with different numbers of scans. The output is an
ordinary Bruker experiment, which you can process in TopSpin or NMRTools.

## From the command line

In the folder that holds the experiments:

```bash
julia -e 'using NMRTools; sumexpts(11, 12, out=999)'
```

adds experiments 11 and 12 and writes the sum as experiment 999. A difference spectrum:

```bash
julia -e 'using NMRTools; sumexpts(11, 12, out=999, weights=[1, -1])'
```

## From Julia

```julia
using NMRTools
cd("myexperiment")

sumexpts(11, 12, out=999)                          # 11 + 12
sumexpts(11, 12, 13, out=999)                      # 11 + 12 + 13
sumexpts(11, 12, out=999, weights=[1, -1])         # 11 − 12
sumexpts(11, 12, 13, out=999, weights=[0.5, 1, 2]) # a weighted sum
```

Experiment numbers refer to folders in the current directory. Paths work too, for experiments
elsewhere:

```julia
sumexpts("data/11", "data/12", out="data/999")
```

## What it does

- The output experiment is a copy of the first input, so it has the same acquisition parameters.
  If the output experiment already exists, `sumexpts` asks before replacing it.
- The raw data (`fid` or `ser`) of the inputs are multiplied by their weights and added.
- The number of scans (`NS`) is set to the weighted sum of the inputs' scans, and the title records
  which experiments were combined and with what weights.
- The processed data copied from the first input are removed, since they no longer match.
- If the experiments have different sizes, they are truncated to the smallest.

`sumexpts` prints a summary as it goes:

```
=== NMR Experiment Summation ===
Output: 999
Inputs: 11, 12
Weights: 1.0, -1.0
Created output experiment from 11
Data format: Float64 (little-endian)
Detected nD experiment, using ser files
Loading 11/ser
Loading 12/ser
Writing output to 999/ser
Scans from 11: 16 × 1.0 = 16.0
Scans from 12: 16 × 1.0 = 16.0
Updated number of scans to 32
=== Operation completed successfully ===
```

## Processing the result

Process the output experiment like any other, in TopSpin or with NMRTools:

```julia
spec = loadfid("999") |> autoprocess
```
