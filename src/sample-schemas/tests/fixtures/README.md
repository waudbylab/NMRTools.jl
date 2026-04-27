# Test fixtures

JSON samples used by the Python, Julia, JavaScript, and MATLAB test suites.
All fixtures are named `sample_v<FROM_VERSION>_<description>.json` and are
expected to migrate to the current schema version (0.4.0) without errors.

| File | From | Purpose |
| --- | --- | --- |
| `sample_v0.0.2_multi.json` | 0.0.2 | Full-tree rename of a 3-component sample; exercises every `rename_key` wildcard path |
| `sample_v0.2.0_multi.json` | 0.2.0 | 3-component sample with `equiv` units to verify wildcard `map` + `set` on arrays |
| `sample_v0.3.0_multi.json` | 0.3.0 | 3-component sample that only needs the v0.4.0 `type` field added |
| `sample_v0.2.0_empty_components.json` | 0.2.0 | Components array present but empty (wildcard ops must be no-ops) |
| `sample_v0.2.0_no_components.json` | 0.2.0 | No `components` key at all (wildcard ops must be silent no-ops) |
| `sample_v0.4.0_already_current.json` | 0.4.0 | Already at latest; migration is a pure no-op |

Post-migration invariants enforced by every language suite:

* `metadata.schema_version == "0.4.0"`
* `metadata.schema_source` ends with `v0.4.0/schema.json`
* `sample.components[*].type` exists on every component (null unless already set)
* `sample.components[*].molecular_weight` exists on every component (added in v0.3.0)
* `sample.components[*]` never retains the `equiv` unit (mapped to `""` in v0.3.0)
* `sample.components[*]` never retains CamelCase keys from v0.0.2 (all renamed)
