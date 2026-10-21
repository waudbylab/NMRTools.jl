# Save spectra and FIDs

[`savenmr`](@ref NMRTools.NMRIO.savenmr) writes a frequency-domain spectrum as Bruker processed
data, and [`savefid`](@ref NMRTools.NMRIO.savefid) writes time-domain data as a raw Bruker
`fid`/`ser`. Both produce experiments that TopSpin opens directly.

```@example save
using NMRTools
fid = exampledata("1D_1H"; fid=true)
spec = autoprocess(fid)
out = mktempdir()
nothing # hide
```

## Spectra

Give the experiment folder and a processing number. If the folder does not exist, it is created,
with the parameter files copied from the experiment the data were loaded from:

```@example save
expt = joinpath(out, "1")
savenmr(expt, spec; procno=2)
readdir(joinpath(expt, "pdata", "2"))
```

`savenmr` writes the real part (`1r`, `2rr`, `3rrr`) and, for multicomplex 1D and 2D data, the
imaginary parts, as 32-bit integers scaled by `NC_proc`. The processing parameter files, both
current (`proc`, `proc2`, …) and status (`procs`, `proc2s`, …), record the axes, the window function
and the phase correction applied. Reading the result back gives the same spectrum:

```@example save
back = loadnmr(joinpath(expt, "pdata", "2"))
maximum(abs, back - realest(spec)) / maximum(abs, back)
```

## FIDs

`savefid` writes time-domain data as a new experiment, for example after linear prediction or
solvent filtering:

```@example save
savefid(joinpath(out, "2"), solventfilter(fid))
readdir(joinpath(out, "2"))[1:5]
```

`savefid` undoes the steps taken by `loadfid`: it restores the first points and the digital-filter
group delay, re-encodes each indirect dimension with its original quadrature scheme and puts the
dimensions back in acquisition order. `loadfid` on the written experiment returns the same data.
Data are written as `Float64` (`DTYPA=2`, as recorded by Avance NEO spectrometers) by default; pass
`dtype=Int32` for older TopSpin versions.

## Parameter files

Both functions copy the acquisition parameter files (`acqus`, `acqu2s`, pulse programme, lists, …)
from the experiment the data were loaded from (its `:experimentfolder` metadata), updating only the
parameters that describe the written data, such as `TD`. If there is no source experiment, minimal
parameter files are generated from the dimension metadata.

Existing data are not overwritten unless `overwrite=true` is passed.

## Opening the result in TopSpin

Open the experiment folder in TopSpin and select the processing number. The spectrum appears as
processed by NMRTools; reprocessing it in TopSpin overwrites it.

!!! note "Sign convention"
    NMRTools spectra are the complex conjugate of TopSpin's: the real parts agree and the imaginary
    parts have the opposite sign. `savenmr` and `loadnmr(...; allcomponents=true)` convert between
    the two, so imaginary data read from or written to Bruker files follow each program's own
    convention. See [Sign convention](../concepts/bruker.md#Sign-convention).

## Peak lists

Peak lists are saved with [`savepeaks`](@ref); see [Save and load peak lists](peaklists.md).

## Pipeline forms

Both functions have forms for the end of a pipeline:

```julia
loadfid("expt/1") |> autoprocess |> savenmr("expt/1"; procno=2)
```
