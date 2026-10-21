# Combining experiments

[`sumexpts`](@ref NMRTools.NMRIO.sumexpts) adds raw Bruker experiments together, with optional
weights, and writes the result as a new experiment. Use it to add repeated acquisitions, to make
difference spectra, or to combine experiments with different numbers of scans. The output is an
ordinary Bruker experiment that can be processed in NMRTools or TopSpin.

```julia
sumexpts(INPUTS...; out=OUTPUT, weights=[...])
```

- `INPUTS...`: experiment numbers (in the current directory) or paths
- `out`: the experiment number or path to write
- `weights`: one weight per input (default: all 1)

It works on 1D (`fid`) and multidimensional (`ser`) experiments.

## Setting up

The examples on this page use a series of ¹⁹F 1D experiments, numbered 1 to 11, from the example
data. `sumexpts` writes into the folder, so we copy the series to a temporary directory first:

```@example sum
using NMRTools, CairoMakie
src = dirname(exampledata("1D_19F_titration"; path=true)[1])
dir = mktempdir()
for i in 1:11
    cp(joinpath(src, string(i)), joinpath(dir, string(i)))
end
expt(n) = joinpath(dir, string(n))
nothing # hide
```

Experiment numbers refer to folders in the current directory, so the examples below run
`sumexpts` inside `cd(dir) do … end`. At the REPL, `cd("myexperiment")` once does the same.

## Adding experiments

Pass the experiment numbers to add, and the experiment number to write:

```@example sum
cd(dir) do
    sumexpts(1, 2; out=101)
end
```

The output experiment is created by copying the first input, so it has the same acquisition
parameters. The number of scans (`NS`) is updated to the weighted sum of the inputs, the title
records which experiments were combined, and the processed data of the copy are removed, since
they no longer match the raw data.

The result is raw data, so process it with [`loadfid`](@ref NMRTools.NMRIO.loadfid) and
[`autoprocess`](@ref):

```@example sum
process(n) = realest(autoprocess(loadfid(expt(n))))
summed = process(101)
summed[:ns], process(1)[:ns]
```

```@example sum
fig, ax, _ = nmrplot(summed[-124.5 .. -123]; title="experiments 1 + 2")
nmrplot!(ax, process(1)[-124.5 .. -123])
nmrplot!(ax, process(2)[-124.5 .. -123])
save("sum-add.svg", fig); nothing # hide
```

![](sum-add.svg)

Paths work as well as experiment numbers, for experiments in another directory:

```julia
sumexpts("myexperiment/11", "myexperiment/12"; out="myexperiment/999")
```

## Difference spectra

Weights of `1` and `-1` subtract one experiment from another:

```@example sum
cd(() -> sumexpts(11, 1; out=102, weights=[1.0, -1.0]), dir)
fig = nmrplot(process(102)[-124.5 .. -123]; title="experiment 11 − experiment 1")
save("sum-diff.svg", fig); nothing # hide
```

![](sum-diff.svg)

Difference spectra are useful for reaction monitoring, comparing spectra before and after a
change, and removing background signals.

## Weighted combinations

Weights scale each input before adding. The number of scans recorded in the output is the
weighted sum:

```@example sum
cd(() -> sumexpts(1, 2, 3; out=103, weights=[0.5, 1.0, 2.0]), dir)
loadfid(expt(103))[:ns]
```

The weights vector must have the same length as the number of inputs. If the experiments have
different sizes, they are truncated to the smallest.

## From the command line

`sumexpts` can be called from a shell, for use in scripts:

```bash
cd myexperiment
julia -e 'using NMRTools; sumexpts(11, 12, out=999)'
julia -e 'using NMRTools; sumexpts(11, 12, out=998, weights=[1.0, -1.0])'
```
