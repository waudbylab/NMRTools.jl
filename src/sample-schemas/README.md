# NMR Sample Schema

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17427433-blue)](https://doi.org/10.5281/zenodo.17427433)

JSON schema definitions for standardising NMR sample metadata across different applications and workflows.

## Overview

This repository contains the authoritative schema definitions for recording metadata about NMR samples, including sample preparation details, experimental conditions, buffer composition, and measurement parameters. The schemas provide a structured, validated format for storing this information alongside experimental data.

## Purpose

Modern NMR experiments generate rich datasets, but associated information about sample preparation (e.g. buffer conditions, protein constructs, labelling schemes) is often recorded inconsistently. These schemas address this by providing:

- Standardised structure for recording sample metadata
- Validation to ensure completeness and consistency
- Version tracking so datasets remain interpretable as the schema evolves
- Interoperability across different tools and applications

## Ecosystem

Currently, sample management is implemented by two applications:

- [NMR Sample Manager (Topspin)](https://nmr-samples.github.io/topspin) - An integrated sample manager for Topspin (v3 upwards)
- [NMR Sample Manager (online)](https://nmr-samples.github.io/online) - A web-based sample manager

Support for sample information is integrated into [NOMAD](https://github.com/nomad-nmr) (v3.6.3 onwards).

Sample parsing is supported by [NMRTools.jl](https://github.com/waudbygroup/NMRTools.jl), and has been integrated into [NMR TITAN](https://www.nmr-titan.com) (development version).

## Schema Structure

The current schema (v0.4.0) organizes sample metadata into the following structure:

| Level 1 | Level 2 | Level 3 | Type / Values |
|---------|---------|---------|---------------|
| **people** | users [array] | | string |
| | groups [array] | | string |
| **sample** | label | | string |
| | physical_form | | enum: `""`, `"solution"`, `"aligned"`, `"solid"` |
| | components [array] | name | string |
| | | type | enum: `""`, `"small molecule"`, `"protein"`, `"protein (intrinsically disordered)"`, `"peptide"`, `"RNA"`, `"DNA"`, `"lipid"`, `"carbohydrate"`, `"other"` |
| | | molecular_weight | number \| null |
| | | concentration_or_amount | number \| null |
| | | unit | enum: `""`, `"uM"`, `"mM"`, `"M"`, `"mg/mL"`, `"%w/v"`, `"%v/v"`, `"mg"`, `"umol"`, `"nmol"` |
| | | isotopic_labelling | enum: `""`, `"natural abundance"`, `"19F"`, `"15N"`, `"13C"`, `"13C,15N"`, `"2H"`, `"2H,15N"`, `"2H,13C,15N"`, `"Ile-13CH3,15N"`, `"ILV-13CH3,15N"`, `"Met-13CH3,15N"`, `"ILVM-13CH3,15N"`, `"2H,Ile-13CH3"`, `"2H,ILV-13CH3"`, `"2H,Met-13CH3"`, `"2H,ILVM-13CH3"`, `"2H,ILVA-13CH3"`, `"2H,ILVMA-13CH3"`, `"2H,ILVMAT-13CH3"`, `"custom"` |
| | | custom_labelling | string |
| **buffer** | ph | | number \| null |
| | components [array] | name | string |
| | | concentration | number \| null |
| | | unit | enum: `""`, `"uM"`, `"mM"`, `"M"`, `"mg/mL"`, `"%w/v"`, `"%v/v"`, `"%w/w"` |
| | chemical_shift_reference | | enum: `""`, `"none"`, `"DSS"`, `"TMS"`, `"TSP"` |
| | reference_concentration | | number \| null |
| | reference_unit | | enum: `""`, `"uM"`, `"mM"`, `"M"`, `"mg/mL"`, `"%w/v"`, `"%v/v"`, `"%w/w"` |
| | solvent | | enum: `""`, `"10% D2O"`, `"100% D2O"`, `"CDCl3"`, `"DMSO-d6"`, `"Methanol-d4"`, `"Acetone-d6"`, `"Acetonitrile-d3"`, `"Benzene-d6"`, `"THF-d8"`, `"custom"` |
| | custom_solvent | | string |
| **nmr_tube** | diameter_mm | | number \| null |
| | type | | enum: `""`, `"regular"`, `"shigemi"`, `"shaped"`, `"coaxial"`, `"J Young"`, `"zirconia rotor"`, `"silicon nitride rotor"`, `"sapphire rotor"` |
| | sample_volume_uL | | number \| null |
| | sample_mass_mg | | number \| null |
| | rack_id | | string |
| | rotor_serial | | string |
| **reference** | sample_id | | string |
| | labbook_entry | | string |
| **notes** | | | string |
| **metadata** | created_timestamp | | string (date-time format) |
| | modified_timestamp | | string (date-time format) |
| | ejected_timestamp | | string (date-time format) |
| | schema_version | | string |
| | schema_source | | string |

## Schema Versions

Schemas are versioned using semantic versioning and tagged in this repository. Each dataset should record the schema version it was created with, ensuring backwards compatibility as the schema evolves.

The schema is currently in the 0.x phase, indicating that the domain model is still being refined through practical use. Breaking changes may occur during this exploratory period. A stable 1.0 release will follow once the model has been validated across multiple use cases.

## Accessing Schemas

Schemas are organised by version, with each version in its own directory:

```
versions/v0.0.1/schema.json
versions/v0.0.2/schema.json
versions/v0.0.3/schema.json
versions/v0.1.0/schema.json
versions/v0.2.0/schema.json
versions/v0.3.0/schema.json
versions/v0.4.0/schema.json
current/schema.json
```

The `current` directory is a copy of the latest tagged release.

To reference a specific schema version in your application:
```
https://github.com/nmr-samples/schema/blob/main/versions/v0.4.0/schema.json
```

To always use the latest schema:
```
https://github.com/nmr-samples/schema/blob/main/current/schema.json
```

## Patching schema updates

The file `current/patch.json` contains methods to update files to the latest schema version. This is written in a simple json DSL:

```json
[
  {
    "from_version": "0.0.2",
    "operations": [
      {"op": "move", "path": "/Users", "to": "/people/users"}
    ]
  }
]
```

| Op | Fields | Behaviour |
|---|---|---|
| `set` | `path`, `value` | Set value at path. On concrete paths, creates intermediate objects if absent. On wildcard paths, missing intermediates are a silent no-op. |
| `remove` | `path` | Remove key at path. No-op if absent. |
| `rename_key` | `path`, `to` | Rename final key segment. No-op if key absent. Error if `to` already exists. |
| `map` | `path`, `from`, `to` | If value at path equals `from`, replace with `to`. Otherwise no-op. |
| `move` | `path`, `to` | Move value to a new path. Creates intermediates. No-op if absent |

Paths: JSON Pointer with `*` wildcard for array elements. Missing intermediate paths → no-op (except `set` which creates them).




## Applications

This schema is used by:

- [NMR Samples (Topspin)](https://nmr-samples.github.io/topspin) - Topspin-integrated sample manager
- [NMR Samples (online)](https://nmr-samples.github.io/online) - Web-based sample manager

## Tests

Unit tests cover the Python, Julia, and JavaScript conversion scripts and run in GitHub Actions (`.github/workflows/test.yml`). MATLAB tests live alongside them but are run locally – MATLAB is proprietary and CI runners are not generally available. Fixtures shared across all four suites are in `tests/fixtures/`.

```
# Python
python -m pytest tests/python

# JavaScript (Node 18+)
cd tests/js && node --test test_migrate.js

# Julia
julia --project=tests/julia -e 'using Pkg; Pkg.instantiate()'
julia --project=tests/julia tests/julia/runtests.jl

# MATLAB (manual, not in CI)
>> cd tests/matlab
>> runtests('testMigrate')
```


## Changelog

### v0.4.0

**Non-breaking changes:**
- Added `sample.components[].type` field (`null | "" | small molecule | peptide | protein | RNA | DNA | carbohydrate | other`) to classify molecule type
- Added `19F` isotopic labelling option: `19F`
- Expanded `buffer.solvent` enum with common NMR solvents: `5% D2O`, `CD2Cl2`, `CD3CN`, `C6D6`, `D6-acetone`, `D5-pyridine`, `D8-toluene`, `D8-THF`, `D12-cyclohexane`, `D3-TFA`

**Infrastructure:**
- Added unit tests for Python, Julia, and JavaScript conversion scripts (CI via GitHub Actions)
- Added MATLAB tests for local runs
- Fixed wildcard `set` operations materialising spurious empty-dict intermediates when the parent array was absent


### v0.3.0

**Breaking Changes:**
- Added `molecular_weight` field (number or null, in Da) to sample components
- Removed `equiv` from sample component `unit` enum
- Renamed `nmr_tube.diameter` to `nmr_tube.diameter_mm`

### v0.2.0

**Non-breaking changes:**
- Moved schema from waudbygroup to new organisation, nmr-samples/schema
- Add `metadata.schema_source` with link to schema

### v0.1.0

**Breaking Changes:**
- Changed `sample.components[].concentration` to `concentration_or_amount` to better reflect that this field can represent either concentration or absolute amounts
- Changed `nmr_tube.diameter` from enum of specific string values ("1.7 mm", "3 mm", "5 mm") to a numeric field (type: number) with min/max validation (0.1-10 mm)
- Renamed `nmr_tube.samplejet_rack_id` to `nmr_tube.rack_id`
- Removed `nmr_tube.samplejet_rack_position` field

**New Features:**
- Added `sample.physical_form` field to specify physical state (solution, aligned, solid)
- Added `nmr_tube.sample_mass_mg` field for recording sample mass
- Added `nmr_tube.rotor_serial` field for solid-state NMR rotors
- Expanded `nmr_tube.type` enum to include rotor types: "zirconia rotor", "silicon nitride rotor", "sapphire rotor"
- Added units for absolute amounts in `sample.components[].unit`: "mg", "umol", "nmol"
- Enhanced `sample.components[].isotopic_labelling` options with deuteration patterns:
  - "2H,15N", "2H,Ile-δ1-13CH3", "2H,Leu/Val-13CH3", "2H,ILV-13CH3"
  - "2H,Met-13CH3", "2H,ILVM-13CH3", "2H,ILVA-13CH3", "2H,ILVMA-13CH3", "2H,ILVMAT-13CH3"
- Reorganized some isotopic labelling options (e.g., separated "Ile-δ1-13CH3,15N" and "ILV-13CH3,15N")

**Description Updates:**
- Updated `people.groups` description to clarify "Research groups (surnames)"
- Changed `nmr_tube` title from "NMR Tube" to "NMR Tube / Rotor"
- Updated `nmr_tube.diameter` description to include rotors
- Updated `nmr_tube.type` description to mention rotors
- Clarified `sample.components[].name` as "Molecule name"
- Improved `sample.components[].concentration_or_amount` description

**Migration Support:**
- Added migration tools in Python and Julia (see `migration-code/`)
- Added `current/patch.json` with automated migration rules
- See [Patching schema updates](#patching-schema-updates) section for migration details

## Contributing

The schema is being developed through practical use in the Waudby laboratory at UCL School of Pharmacy. If you're using these schemas in your own work and have suggestions for improvements or extensions, please open an issue to discuss or contact [Chris](mailto:c.waudby@ucl.ac.uk) directly. We're especially interested in hearing from groups who might want to adopt or extend these schemas for their own workflows.
