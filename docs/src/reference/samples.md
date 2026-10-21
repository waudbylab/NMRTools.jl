# Sample metadata schema

Sample metadata are JSON files describing the physical sample, following the
[NMR sample schema](https://nmrsamples.github.io). NMRTools reads the version bundled with it
(v0.4.0) and migrates files written with earlier versions. A sample file is associated with an
experiment when the experiment's date falls between the file's `created_timestamp` and
`ejected_timestamp`; sample files are looked for in the folder that contains the experiment folders.

Fields are read with [`sample`](@ref), passing the section and field names; see
[Find sample information](../howto/samples.md).

When a field such as `solvent` or `isotopic_labelling` has the value `"custom"`, `sample` returns the
corresponding `custom_…` field instead.

## `people`

People

| Field | Type | Description |
|:--|:--|:--|
| `users` | array | People involved in this experiment |
| `groups` | array | Research groups (surnames) involved in this experiment |

## `sample`

Sample

| Field | Type | Description |
|:--|:--|:--|
| `label` | string | Quick sample identifier |
| `physical_form` | string | Physical form of sample |
| `components` | array | Sample components (allows for multiple components) |

## `buffer`

Buffer

| Field | Type | Description |
|:--|:--|:--|
| `ph` | number | Buffer pH |
| `components` | array | Buffer components and their concentrations |
| `chemical_shift_reference` | string | Chemical shift reference compound |
| `reference_concentration` | number | Concentration of chemical shift reference |
| `reference_unit` | string | Unit for reference concentration |
| `solvent` | string | Solvent composition |
| `custom_solvent` | string | Custom solvent description if 'custom' is selected above |

## `nmr_tube`

NMR tube or rotor properties

| Field | Type | Description |
|:--|:--|:--|
| `diameter_mm` | number | NMR tube or rotor diameter (mm) |
| `type` | string | Type of NMR tube or rotor |
| `sample_volume_uL` | number | Sample volume (µL) |
| `sample_mass_mg` | number | Sample mass (mg) |
| `rack_id` | string | Identifier for SampleJet rack |
| `rotor_serial` | string | Rotor serial number |

## `reference`

Cross-references for lab notebooks etc.

| Field | Type | Description |
|:--|:--|:--|
| `sample_id` | string | Local sample identifier |
| `labbook_entry` | string | Reference to laboratory notebook entry |

## `notes`

Free text field for additional notes and observations

Type: string

## `metadata`

Metadata

| Field | Type | Description |
|:--|:--|:--|
| `created_timestamp` | string | ISO 8601 timestamp when this record was created |
| `modified_timestamp` | string | ISO 8601 timestamp when this record was last modified |
| `ejected_timestamp` | string | ISO 8601 timestamp when this sample was ejected |
| `schema_version` | string | Version of this schema used |
| `schema_source` | string | Source of this schema version |

