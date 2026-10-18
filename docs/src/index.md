# NMRTools.jl

NMRTools.jl is a Julia library for NMR spectroscopy data. It reads and writes 1D, 2D and higher-dimensional datasets, processes raw time-domain data, and exposes spectra as named-dimension arrays with full metadata support. Quadrature is carried by multicomplex numbers, so each dimension of a multidimensional experiment has its own imaginary unit.

- **Multi-format support**: read Bruker (processed and raw `fid`/`ser`), nmrPipe and UCSF/Sparky data; write Bruker processed data and raw FIDs
- **Processing**: composable apodization, zero-filling, Fourier transform, phasing, linear prediction and solvent suppression, with automatic phasing and a one-step `autoprocess`
- **Peaks**: N-dimensional peak picking, peak lists that are Tables.jl tables, registration of a peak list onto a spectrum, and NMRAnalysis.jl and Sparky peak list files
- **Assignments**: import assigned chemical shifts from the BMRB and build peak lists from them
- **Chemical shift indexing**: array-like indexing with chemical shift values using `spec[8.0 .. 9.0]` syntax
- **Plotting**: publication-ready plots via Plots.jl recipes or Makie, with contour level and colour control
- **Metadata**: acquisition parameters from Bruker `acqus` files, pulse programme annotations and sample information; window functions, power levels and frequency lists represented as typed objects
- **DimensionalData.jl integration**: named dimension arrays with selectors for frequency, time and gradient dimensions

## Installation

Install NMRTools.jl through the Julia package manager:

```julia
using Pkg
Pkg.add("NMRTools")
```

Load and plot a spectrum:

```julia
using NMRTools, Plots
spec = loadnmr("path/to/experiment/1")
plot(spec)
```

See the [Getting started](@ref) guide for more examples, the tutorials for loading, processing, plotting and peak picking, and the worked examples for complete relaxation and diffusion analyses.


## Related packages

NMRTools.jl is part of a suite of Julia packages for NMR data handling developed
by the [Waudby lab](https://waudbylab.org):

- **[MulticomplexNumbers.jl](https://waudbylab.org/MulticomplexNumbers.jl)** — Multicomplex number types with FFTW extension, providing the natural algebraic representation of multi-dimensional NMR data.

- **[NMRAnalysis.jl](https://waudbylab.org/NMRAnalysis.jl)** — Quantitative NMR data analysis, with a particular focus on biomolecular relaxation, diffusion, exchange, and protein dynamics measurements.

- **[NMRScreen.jl](https://waudbylab.org/NMRScreen.jl)** — tools for
  fragment screening by NMR, including automated analysis of large compound libraries.
