# Pulse programme annotations

!!! warning "Alpha Development Status"
    This annotation system is in early alpha development (currently v0.0.2). The schema, syntax, and functionality are subject to significant changes. Use with caution in production environments.

## Motivation

Modern NMR pulse programmes are cryptic text files that encode complex experimental procedures, but they lack semantic information about what they actually measure. Analysis software must guess at the meaning of parameters like `VALIST` or `VPLIST`, leading to errors and requiring manual interpretation.

This annotation system embeds structured semantic metadata directly within pulse programme files as comments, enabling:

- **Automated Analysis**: Software can automatically understand experimental structure
- **Reproducible Research**: Clear provenance and citation information
- **Knowledge Preservation**: Experimental design captured alongside implementation
- **Cross-Platform Compatibility**: Vendor-agnostic semantic layer
- **Collaborative Development**: Track contributors and development status

!!! note "Current Implementation"
    Basic annotation parsing is implemented in NMRTools.jl. The system is under active development with ongoing refinement of the schema and parsing capabilities.

## Design Principles

### Human and Machine Readable
Annotations use simple comment syntax that NMR spectroscopists can read naturally whilst being structured enough for automated parsing.

### Embedded and Immutable
Metadata travels with the pulse programme file and cannot drift out of sync with the code.

### Modular and Extensible
Experiments are described as combinations of core types and features, allowing new variations without redefining the entire taxonomy.

### Version Controlled
Both individual pulse programmes and the annotation system itself are versioned for reproducibility and compatibility.

## Annotation Syntax

### Basic Format

Annotations are extensions of Bruker comment lines, marked with `;@ `.
Within these lines, annotations are written in YAML format. The initial `;@ ` will be
stripped when parsing, permitting multi-line entries.
```
;@ parameter: text value
;@ parameter1: 0.123
;@ parameter2: [list, of, items]
;@ parameter3: {associative: array, key: value}
;@ parameter4: |
;@   This is
;@   a multi-line
;@   entry.
```

### Parameter References

Parameter and list variable names are referenced directly without special syntax:
```yaml
;@ r1rho:
;@   power: VALIST      # List variable
;@   duration: VPLIST   # List variable  
;@   offset: cnst28     # Constant parameter
```

The context distinguishes parameter names (letters) from literal values (numbers). NMRTools automatically resolves these references to actual values from the acqus file when loading data.

### Dimension References

Dimensions use dotted path notation that references into experiment-specific parameter blocks:
```yaml
;@ dimensions: [r1rho.power, r1rho.duration, f1]
;@ r1rho:
;@   power: VALIST
;@   duration: VPLIST
```

This creates an explicit connection: `r1rho.power` refers to the `power` field in the `r1rho` block, which is controlled by the `VALIST` parameter.

## Using Annotations in NMRTools.jl

The complete annotation schema is defined in the [pulse programmes repository](https://github.com/waudbylab/pulseprograms). See the [schema documentation](https://waudbylab.org/pulseprograms/schema/fields/) and [controlled vocabulary](https://github.com/waudbylab/pulseprograms/blob/main/VOCABULARY.md) for details.

In NMRTools.jl, parsed annotation data is accessible via the `:annotations` metadata field or the [`annotations()`](@ref) convenience function.

### Accessing annotation data

Annotations can be accessed using the `annotations()` function, which provides convenient nested access to annotation data. Dotted notation is automatically expanded:
```julia
# Load a spectrum with annotations
spec = loadnmr("path/to/annotated/experiment")

# Access all annotations
all_annotations = annotations(spec)

# Access specific annotation fields
experiment_type = annotations(spec, "experiment_type")

# Access nested dictionary fields - these are equivalent:
r1rho_duration = annotations(spec, "r1rho", "duration")
r1rho_duration = annotations(spec, "r1rho.duration")

# Access array elements by index
first_dimension = annotations(spec, "dimensions", 1)

# Deep nesting with dotted notation
calibration_start = annotations(spec, "calibration.duration.start")
```

The `annotations()` function accepts both string and symbol keys, and returns `nothing` if the requested field does not exist.

### Helper Functions

#### `referencepulse(spec, nucleus)`

Get the reference pulse calibration for a given nucleus:
```julia
# Get 19F reference pulse (returns tuple of pulse length and power)
p1, pl1 = referencepulse(spec, "19F")

# Get 1H reference pulse
p3, pl2 = referencepulse(spec, :1H)
```

Returns `nothing` if no reference pulse is found for the specified nucleus.

### Example: 19F R1ρ on-resonance experiment

This example experiment is annotated as follows (schema v0.0.2):
```
;@ schema_version: "0.0.2"
;@ sequence_version: "0.1.0"
;@ title: 19F on-resonance R1rho relaxation dispersion
;@ authors:
;@   - Chris Waudby <c.waudby@ucl.ac.uk>
;@   - Jan Overbeck
;@ created: 2020-01-01
;@ last_modified: 2025-11-15
;@ repository: github.com/waudbylab/pulseprograms
;@ status: beta
;@ experiment_type: [r1rho, 1d]
;@ features: [on_resonance, temperature_compensation]
;@ typical_nuclei: [19F, 1H]
;@ citation:
;@   - Overbeck (2020)
;@ dimensions: [r1rho.power, r1rho.duration, f1]
;@ acquisition_order: [f1, r1rho.duration, r1rho.power]
;@ reference_pulse:
;@   - {channel: f1, pulse: p1, power: pl1}
;@   - {channel: f2, pulse: p3, power: pl2}
;@ r1rho:
;@   channel: f1
;@   power: VALIST
;@   duration: VPLIST
;@   offset: 0
;@   alignment: hard_pulse
```

Using annotated data, we can access the experiment metadata:
```julia
using NMRTools

# Load the annotated experiment
spec = loadnmr("path/to/experiment")

# Access experiment metadata
annotations(spec, "title")           # "19F on-resonance R1rho relaxation dispersion"
annotations(spec, "experiment_type") # ["r1rho", "1d"]
annotations(spec, "features")        # ["on_resonance", "temperature_compensation"]

# Access nested r1rho parameters (dotted notation works)
annotations(spec, "r1rho.power")     # Vector of power values resolved from VALIST
annotations(spec, "r1rho.duration")  # Vector of duration values resolved from VPLIST
annotations(spec, "r1rho", "channel") # "f1"
annotations(spec, "r1rho", "offset")  # 0 (on-resonance)

# Access dimension information
annotations(spec, "dimensions")       # ["r1rho.power", "r1rho.duration", "f1"]
annotations(spec, "dimensions", 1)    # "r1rho.power"

# Get reference pulse calibration
p1, pl1 = referencepulse(spec, "19F")
```

### Additional Examples

#### CEST Experiment
```yaml
;@ schema_version: "0.0.2"
;@ sequence_version: "1.0.0"
;@ title: 19F CEST
;@ experiment_type: [cest, 1d]
;@ typical_nuclei: [19F]
;@ dimensions: [cest.offset, f1]
;@ acquisition_order: [f1, cest.offset]
;@ reference_pulse:
;@   - {channel: f1, pulse: p1, power: pl1}
;@ cest:
;@   channel: f1
;@   power: pl8
;@   duration: d18
;@   offset: FQ1LIST
```

#### Diffusion Experiment
```yaml
;@ schema_version: "0.0.2"
;@ sequence_version: "2.0.1"
;@ title: 1H STE diffusion
;@ experiment_type: [diffusion, 1d]
;@ features: [ste, watergate]
;@ typical_nuclei: [1H]
;@ dimensions: [diffusion.gradient_strength, f1]
;@ acquisition_order: [f1, diffusion.gradient_strength]
;@ reference_pulse:
;@   - {channel: f1, pulse: p1, power: pl1}
;@ diffusion:
;@   type: bipolar
;@   coherence: [f1, 1]
;@   big_delta: d20
;@   little_delta: p31
;@   tau: d17
;@   gradient_strength: {type: linear, start: cnst1, end: cnst2, scale: gpz6}
;@   gradient_shape: gpnam6
```

## Resources

- [Pulse Programme Repository](https://github.com/waudbylab/pulseprograms)
- [Schema Documentation](https://waudbylab.org/pulseprograms/schema/fields/)
- [Controlled Vocabulary (VOCABULARY.md)](https://github.com/waudbylab/pulseprograms/blob/main/VOCABULARY.md)
- [Decision Log](https://github.com/waudbylab/pulseprograms/blob/main/DECISIONS.md) - Design rationale and schema evolution

## Automatic Dimension Annotation

When loading annotated data with `loadnmr`, NMRTools automatically applies semantic dimension types based on the `dimensions` array in the annotations. This eliminates manual calls to dimension setter functions.

### Supported Dimension Types

| Annotation Pattern | Dimension Type | Units | Example |
|-------------------|----------------|-------|---------|
| `*.duration` | `TrelaxDim` | seconds | `r1rho.duration`, `relaxation.duration` |
| `*.power` | `SpinlockDim` | Hz | `r1rho.power` |
| `*.offset` | `OffsetDim` | ppm | `cest.offset` |
| `*.g` | `G1Dim`/`G2Dim` | T/m | `diffusion.g` |

### Automatic Conversions

The annotation system automatically performs unit conversions:

- **Power → Hz**: Power values (from `VALIST`) are converted to RF field strength in Hz using the reference pulse calibration
- **FQList → ppm**: Frequency lists (from `FQ1LIST` etc.) are converted from Hz offsets to absolute ppm values. The `channel` of the block is resolved against the spectrometer [channel model](@ref "Spectrometer channels"): if the channel's nucleus has a detected axis, its referencing is used; otherwise the channel's own base/carrier frequencies are used. This means an offset list still converts to ppm when the manipulated nucleus is never detected — for example a ¹³C CEST profile read out on a ¹H axis.
- **Gradient scaling**: Gradient values are scaled by `gmax` to give actual field strengths

### Example

```julia
# Load annotated R1rho data
spec = loadnmr("path/to/r1rho/experiment")

# Dimensions are automatically typed based on annotations
dims(spec)
# (F1Dim(...), SpinlockDim(...), TrelaxDim(...))

# Access dimension values directly
collect(dims(spec, SpinlockDim))  # RF field strengths in Hz
collect(dims(spec, TrelaxDim))   # Spinlock durations in seconds

# Metadata is also set
spec[SpinlockDim, :units]  # "Hz"
spec[TrelaxDim, :units]    # "s"
spec[TrelaxDim, :delay_type]  # :r1rho
```

### Manual Dimension Setting

If annotations are not available, dimensions can be set manually using:

- `setrelaxtimes` - Set relaxation/delay dimension
- `setspinlockfield` - Set spinlock field strength dimension
- `setoffsets` - Set offset dimension
- `setgradientlist` - Set gradient dimension