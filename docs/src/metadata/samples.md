# Sample metadata

Sample metadata can be automatically loaded when present, following the schema defined at [nmr-sample-schema](https://github.com/waudbylab/nmr-sample-schema). This provides structured information about the physical sample, buffer composition, NMR tube characteristics, and personnel involved in the experiment.

Sample metadata is stored as JSON files in the experiment directory and is automatically associated with spectra based on timestamps:

```julia
# Load spectrum — sample metadata loaded automatically if present
spec = loadnmr("path/to/experiment")

# Check if sample metadata is available
hassample(spec)

# Get the path to the matched sample JSON file
samplefile(spec)
```

## Accessing sample metadata

Use the [`sample`](@ref) function to access sample metadata:

```julia
# Get the NMRSample object
sample(spec)
```

This returns an [`NMRSample`](@ref) wrapping the complete sample metadata structure. Navigate nested fields by passing keys:

```julia
# Get sample label
sample(spec, :sample, :label)

# Get user names
sample(spec, "people", "users")

# Get buffer solvent
sample(spec, :buffer, :solvent)

# Get component list
sample(spec, "sample", "components")

# Get component names (key mapped over array elements)
sample(spec, :sample, :components, :name)
```

Keys can be strings or symbols and are case-insensitive. If a key is not found at any level, `nothing` is returned.

When a key is applied to an array (such as `components`), it is mapped over each element and a filtered array of results is returned.

For `isotopic_labelling` and `solvent`, if the stored value is `"custom"`, the corresponding `custom_labelling` or `custom_solvent` field is returned automatically:

```julia
# Returns "custom_labelling" value if isotopic_labelling == "custom"
sample(spec, :sample, :components, :isotopic_labelling)

# Returns "custom_solvent" value if solvent == "custom"
sample(spec, :buffer, :solvent)
```

## Schema structure

The sample metadata follows a hierarchical structure with the following main sections:

- **sample**: Sample identification and component list
  - `label`: Descriptive label for the sample
  - `components`: Array of sample components
    - `name`, `isotopic_labelling`, `concentration`, `unit`
- **buffer**: Buffer composition
  - `solvent`: Solvent description (e.g., "10% D2O")
  - `chemical_shift_reference`: Reference compound used
  - `components`: Array of buffer components
    - `name`, `concentration`, `unit`
- **nmr_tube**: NMR tube specifications
  - `type`: Tube type (e.g., "regular", "Shigemi")
  - `diameter`: Tube diameter (e.g., "5 mm")
  - `sample_volume_uL`: Sample volume in microliters
- **people**: Personnel information
  - `users`: Array of user names
  - `groups`: Array of group names
- **metadata**: Schema and timestamp information
  - `schema_version`: Version of the schema used
  - `created_timestamp`: When the sample was inserted
  - `ejected_timestamp`: When the sample was removed from the spectrometer
  - `modified_timestamp`: When the sample metadata was last modified
- **notes**: Free-text notes about the sample

For the complete schema specification, see the [nmr-sample-schema repository](https://github.com/waudbylab/nmr-sample-schema).

## Timestamp matching

Sample metadata files are matched to experiments based on timestamps:

1. The experiment acquisition date is read from the `acqus` file
2. Sample JSON files in the parent directory are scanned
3. A sample is matched if the acquisition date falls between the sample's `created_timestamp` and `ejected_timestamp`
4. If no `ejected_timestamp` is present, the sample is assumed to still be in the spectrometer

## Scanning directories

For workflows that process many experiments at once, use the scanning API to avoid loading full binary data:

```julia
# Scan all experiments in a directory (metadata only, no binary data loaded)
expts = scanexperiments("/nmr/projects/lysozyme/")

# Check sample association
hassample(expts[1])
sample(expts[1], "sample", "label")

# Scan all sample files in the same directory (sample JSONs live alongside experiment folders)
samples = scansamples("/nmr/projects/lysozyme/")

# Find the sample matching a specific experiment
s = findsample(expts[1])
s = findsample(expts[1], samples)   # faster — uses pre-scanned list

# Find all experiments for a given sample
group = findexperiments(samples[1], expts)
```
