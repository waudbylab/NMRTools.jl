# Find acquisition parameters

When a spectrum or FID is loaded, the Bruker acquisition files (`acqus`, `acqu2s`, …), the pulse
programme and any lists in the experiment folder are read. [`acqus`](@ref) looks up their values.

```@setup acq
using NMRTools
spec = exampledata("2D_HN")
relax = exampledata("pseudo2D_T2")
```

## Parameters

Give `acqus` a Bruker parameter name, as a string in any case or as a lowercase symbol. Arrayed
parameters (pulse lengths, delays, powers, …) take an index, numbered as in TopSpin. Times are
converted to seconds and frequencies to Hz (see [Time and frequency units](../concepts/units.md)):

```@repl acq
spec = exampledata("2D_HN");
acqus(spec, :pulprog)
acqus(spec, "TE")          # temperature, K
acqus(spec, :ns)
acqus(spec, :p, 1)         # pulse length P1, s
acqus(spec, :d, 1)         # delay D1, s
```

Without an index, an arrayed parameter is returned whole, as a dictionary keyed by the TopSpin
index (so `acqus(spec, :d)[0]` is `D0`), which avoids confusing Bruker's zero-based arrays with
Julia's one-based ones. The parameters of indirect dimensions, from `acqu2s`, `acqu3s`, …, are in
the spectrum metadata under `:acqu2s`, `:acqu3s`, …, and the values that describe each dimension
are also in its own metadata:

```@repl acq
spec[:acqu2s][:td]
spec[F2Dim, :td], spec[F2Dim, :swhz]
```

## Powers and channels

Power levels are returned as [`Power`](@ref) objects, which convert between dB and W (see
[Power levels](../reference/power.md), which also covers RF field strengths). The RF channels
(`:f1`…`:f8`) are described by `channels`, which gives the frequencies of every nucleus pulsed in
the experiment, including those without a detected axis, such as the carbon channel of an
experiment read out on ¹H and ¹⁵N:

```@repl acq
p = acqus(spec, :plw, 1);
watts(p), db(p)
channels(spec, :f1)[:nucleus]
channels(spec, "13C")[:bf]     # Hz
```

A channel can be looked up by label (`:f1`), by `Nucleus` (`C13`) or by a string (`"13C"`), and its
dictionary has keys `:label`, `:nucleus`, `:bf` (base frequency), `:sfo` (carrier frequency),
`:offsethz` and `:offsetppm`.

## Lists and other files

`vclist`, `vdlist`, `valist`, `vplist` and frequency lists (`fq1list`…) are read with `acqus` too,
in NMRTools units: delays and pulse lengths in seconds, amplitudes in dB, and frequency lists as
[`FQList`](@ref) objects (see [Frequency lists](../reference/fqlists.md)). `acqus` returns
`nothing` for a list that is not present.

```@repl acq
relax = exampledata("pseudo2D_T2");
acqus(relax, :vclist)
```

Pulse programmes written with NMRTools' annotation format describe their own parameters; see
[Pulse programme annotations](../experimental/annotations.md).
[`loadmetadata`](@ref NMRTools.NMRIO.loadmetadata) reads all of this from an experiment folder
without loading the data, and [`scanexperiments`](@ref NMRTools.NMRIO.scanexperiments) does so for
every experiment in a folder (see [Find sample information](samples.md)).
