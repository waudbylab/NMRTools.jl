| Op | Fields | Behaviour |
|---|---|---|
| `set` | `path`, `value` | Set value at path. On concrete paths, creates intermediate objects if absent. On wildcard paths, missing intermediates are a silent no-op (no empty containers are materialised). |
| `remove` | `path` | Remove key at path. No-op if absent. |
| `rename_key` | `path`, `to` | Rename final key segment. No-op if key absent. Error if `to` already exists. |
| `map` | `path`, `from`, `to` | If value at path equals `from`, replace with `to`. Otherwise no-op. |
| `move` | `path`, `to` | Move value to a new path. Creates intermediates. No-op if absent |

Paths: JSON Pointer with `*` wildcard for array elements. Missing intermediate paths → no-op (except `set` which creates them).
