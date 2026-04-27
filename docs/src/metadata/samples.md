# Sample metadata

Sample metadata can be automatically loaded when present, following the schema defined at [nmrsamples.github.io](https://nmrsamples.github.io). This provides structured information about the physical sample, buffer composition, NMR tube characteristics, and personnel involved in the experiment.

Sample metadata is stored as JSON files in the experiment directory and is automatically associated with spectra based on timestamps:

```julia
# Load spectrum — sample metadata loaded automatically if present
spec = loadnmr("path/to/experiment")

# Check if sample metadata is available
hassample(spec)
# true

# Get the path to the matched sample JSON file
samplefile(spec)
# "/nmr/projects/lysozyme/2025-08-23_163924_lysozyme.json"
```

## Accessing sample metadata

Use the [`sample`](@ref) function to access sample metadata:

```julia
# Get the NMRSample object
sample(spec)
# NMRSample("/nmr/projects/lysozyme/2025-08-23_163924_lysozyme.json")
```

Navigate nested fields by passing keys:

```julia
# Get sample label
sample(spec, :sample, :label)
# "lysozyme"

# Get user names
sample(spec, "people", "users")
# ["jsmith"]

# Get buffer solvent
sample(spec, :buffer, :solvent)
# "95% H2O / 5% D2O"

# Get component list
sample(spec, "sample", "components")
# [Dict("name" => "HEWL", "concentration" => 0.69, "unit" => "mM", ...)]

# Get component names (key mapped over array elements)
sample(spec, :sample, :components, :name)
# ["HEWL", "gadodiamide"]
```

Keys can be strings or symbols and are case-insensitive. If a key is not found at any level, `nothing` is returned.

When a key is applied to an array (such as `components`), it is mapped over each element and a filtered array of results is returned.

For `isotopic_labelling` and `solvent`, if the stored value is `"custom"`, the corresponding `custom_labelling` or `custom_solvent` field is returned automatically:

```julia
# Returns "custom_labelling" value if isotopic_labelling == "custom"
sample(spec, :sample, :components, :isotopic_labelling)
# ["2H,13C-Met", "unlabelled"]

# Returns "custom_solvent" value if solvent == "custom"
sample(spec, :buffer, :solvent)
# "95% H2O / 5% D2O + glycerol"
```

## Schema

For the complete schema specification, including all available fields and their meanings, see [nmrsamples.github.io](https://nmrsamples.github.io).

## Scanning directories

For workflows that process many experiments at once, use the scanning API to avoid loading full binary data:

```julia
# Scan all experiments in a directory (metadata only, no binary data loaded)
expts = scanexperiments("/nmr/projects/lysozyme/")

# Check sample association
hassample(expts[1])
# true

sample(expts[1], "sample", "label")
# "lysozyme"

# Scan all sample files in the same directory (sample JSONs live alongside experiment folders)
samples = scansamples("/nmr/projects/lysozyme/")

# Find the sample matching a specific experiment
s = findsample(expts[1])
s = findsample(expts[1], samples)   # faster — uses pre-scanned list

# Find all experiments for a given sample
group = findexperiments(samples[1], expts)
```
