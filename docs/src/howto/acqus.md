# Find acquisition parameters

When a spectrum or FID is loaded, the Bruker acquisition files (`acqus`, `acqu2s`, …), the pulse
programme and any lists in the experiment folder are read. [`acqus`](@ref) looks up their values.

```@example acq
using NMRTools
spec = exampledata("2D_HN")
nothing # hide
```

## Parameters

Give a parameter name, as a string in any case (`"TE"`, `"te"`) or a lowercase symbol (`:te`):

```@example acq
acqus(spec, :te), acqus(spec, "PULPROG"), acqus(spec, :ns)
```

Times are converted to seconds and frequencies to Hz (see
[Time and frequency units](../concepts/units.md)).

## Arrayed parameters

Pulse lengths, delays, powers and other arrayed parameters take an index, as numbered in TopSpin:

```@example acq
acqus(spec, :p, 1), acqus(spec, :d, 1)
```

Without an index the whole array is returned, as a dictionary keyed by the TopSpin index, so `d[0]`
is `D0`:

```@example acq
d = acqus(spec, :d)
d[0], d[1]
```

Arrays are dictionaries rather than vectors to avoid confusion between Bruker arrays, which count
from zero, and Julia arrays, which count from one.

## Indirect dimensions

Parameters from `acqu2s`, `acqu3s`, … are under the keys `:acqu2s`, `:acqu3s` of the spectrum
metadata:

```@example acq
spec[:acqu2s][:td]
```

Most dimension-specific values are also in each dimension's metadata, for example
`spec[F2Dim, :td]` or `spec[F2Dim, :swhz]`.

## Power levels

Power levels are returned as [`Power`](@ref) objects, which convert between dB and W:

```@example acq
p = acqus(spec, :plw, 1)
watts(p), db(p)
```

See [Power levels](../reference/power.md), including conversion to RF field strength.

## Spectrometer channels

Frequency dimensions carry their own `:bf`, `:sf` and `:nucleus` metadata, but only for detected
nuclei. To reach a nucleus that was pulsed but never detected, such as the carbon channel of a
¹³C-edited experiment read out on ¹H, use the channel model. `channels(spec)` returns all the RF
channels (`:f1`…`:f8`), and `channels(spec, ...)` one channel, by label, by `Nucleus` or by a
string:

```@example acq
channels(spec, :f1)
```

```@example acq
channels(spec, "13C")[:bf]     # carbon base frequency, with no carbon axis
```

Each channel dictionary has keys `:label`, `:nucleus`, `:bf` (base frequency, Hz), `:sfo` (carrier
frequency, Hz), `:offsethz` and `:offsetppm`.

## Lists

`vclist`, `vdlist`, `valist`, `vplist` and frequency lists (`fq1list`…) in the experiment folder are
read too:

```@example acq
relax = exampledata("pseudo2D_T2")
acqus(relax, :vclist)
```

Values are converted to NMRTools units: delays and pulse lengths in seconds, amplitudes in dB, and
frequency lists as [`FQList`](@ref) objects (see [Frequency lists](../reference/fqlists.md)).
`acqus` returns `nothing` for a list that is not present.

## Pulse programme annotations

Pulse programmes annotated with NMRTools' annotation format describe their own parameters; see
[Pulse programme annotations](../experimental/annotations.md).

## Without loading the data

[`loadmetadata`](@ref NMRTools.NMRIO.loadmetadata) reads the acquisition metadata of an experiment
folder without loading any spectrum, and [`scanexperiments`](@ref NMRTools.NMRIO.scanexperiments)
does this for every experiment in a folder (see [Find sample information](samples.md)).
