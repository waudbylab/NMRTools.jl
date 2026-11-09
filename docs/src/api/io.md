# Loading and saving

```@index
Pages = ["io.md"]
```

## Loading

```@docs
loadnmr
loadfid
loadmetadata
loadjdx
```

## Example data

NMRTools re-exports these from [NMRExampleData.jl](https://waudbylab.org/NMRExampleData.jl/stable/),
whose [API reference](https://waudbylab.org/NMRExampleData.jl/stable/api/) documents them:

- `example(name)`: the path or paths to load, downloading the dataset the first time;
  `example(name, format)` for the same data in another format, such as `:nmrpipe` or `:ucsf`
- `examples()`: the names of the datasets; `examples(format)` those with a given format
- `exampleinfo(name)`: a dataset's title, description, credit, licence and formats
- `exampleanalysis(name)`: the NMRAnalysis analysis saved with a dataset

## Saving

```@docs
savenmr
savefid
```

## Combining experiments

```@docs
sumexpts
```
