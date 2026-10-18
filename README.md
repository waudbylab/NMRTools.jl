# NMRTools.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://waudbylab.org/NMRTools.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://waudbylab.org/NMRTools.jl/dev)
[![CI](https://github.com/waudbylab/NMRTools.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/waudbylab/NMRTools.jl/actions/workflows/Runtests.yml)
[![Codecov](https://codecov.io/gh/waudbylab/NMRTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/waudbylab/NMRTools.jl)
[![DOI](https://zenodo.org/badge/251587402.svg)](https://zenodo.org/badge/latestdoi/251587402)

![NMRTools logo](logo.png)

NMRTools.jl is a Julia library for NMR spectroscopy data. It reads and writes 1D, 2D and higher-dimensional datasets, processes raw time-domain data, and exposes spectra as named-dimension arrays with full metadata support. Quadrature is carried by multicomplex numbers, so each dimension of a multidimensional experiment has its own imaginary unit.

## Features

- **Multi-format support**: read Bruker (processed and raw `fid`/`ser`), nmrPipe and UCSF/Sparky data; write Bruker processed data and raw FIDs
- **Processing**: composable apodization, zero-filling, Fourier transform, phasing, linear prediction and solvent suppression, with automatic phasing and a one-step `autoprocess`
- **Peaks**: N-dimensional peak picking, peak lists that are Tables.jl tables, registration of a peak list onto a spectrum, and NMRAnalysis.jl and Sparky peak list files
- **Assignments**: import assigned chemical shifts from the BMRB and build peak lists from them
- **Chemical shift indexing**: array-like indexing with chemical shift values using `spec[8.0 .. 9.0]` syntax
- **Plotting**: publication-ready plots via Plots.jl recipes or Makie, with contour level and colour control
- **Metadata**: acquisition parameters from Bruker `acqus` files, pulse programme annotations and sample information; window functions, power levels and frequency lists represented as typed objects
- **DimensionalData.jl integration**: named dimension arrays with selectors for frequency, time and gradient dimensions

## Quick Start

Install NMRTools.jl through the Julia package manager:

```julia
using Pkg
Pkg.add("NMRTools")
```

Load and plot a processed spectrum:

```julia
using NMRTools, Plots

# Load data (auto-detects format)
spec = loadnmr("path/to/experiment")

# Plot the full spectrum
plot(spec)

# Zoom to a chemical shift range
plot(spec[8.0 .. 9.0])

# Access data and metadata
data(spec)
metadata(spec)
```

Process raw data, either in one step or with explicit control of each stage:

```julia
spec = loadfid("path/to/experiment") |> autoprocess

spec = loadfid("path/to/hsqc") |> apodize(Cos2Window()) |> zerofill |> ft |>
       phase(F1Dim, -50) |> phase(F2Dim, 0) |> real
```

Pick peaks and save them for analysis in NMRAnalysis.jl:

```julia
pks = detectpeaks(spec; snr=20)
savepeaks("peaklist.csv", pks)
```

## Documentation

Documentation with tutorials and examples is available at:
- [**Stable docs**](https://waudbylab.org/NMRTools.jl/stable) — Latest release
- [**Dev docs**](https://waudbylab.org/NMRTools.jl/dev) — Development version

## Development Status

> **NOTE**: This package is under active development and the API may change between minor releases. Loading and processing raw time-domain data is new in v0.3 and its interface is the most likely to change. Please refer to the [documentation](https://waudbylab.org/NMRTools.jl/stable) for current usage examples.

## Citation

If you use NMRTools.jl in your research, please cite:
[![DOI](https://zenodo.org/badge/251587402.svg)](https://zenodo.org/badge/latestdoi/251587402)
