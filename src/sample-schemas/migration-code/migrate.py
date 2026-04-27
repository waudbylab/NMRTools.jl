"""
update_to_latest_schema(data, migrations_path=None)

Migrate NMR sample metadata to the latest schema version.

Call update_to_latest_schema(data) with a parsed JSON dict. It modifies
the dict in place and returns it.

The migration patch file is expected at ../current/patch.json relative to
this module. Override by setting _MIGRATIONS_PATH or passing migrations_path.

Compatible with Jython 2.7 and CPython 2.7+/3.x.
"""

import json
import os

# set this to the path of the patch file (by default, looks in ../current/patch.json relative to this script)
_MIGRATIONS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "current", "patch.json")


def _parse_path(path):
    """Split a JSON Pointer string into segments, handling '/' and '~' escaping."""
    if not path:
        return []
    if not path.startswith("/"):
        raise ValueError("Path must start with '/': " + path)
    parts = path[1:].split("/")
    return [p.replace("~1", "/").replace("~0", "~") for p in parts]


def _resolve(data, segments):
    """Return list of (parent, key) pairs for all locations matching a path with wildcards.

    Each pair satisfies: parent[key] is the value at that location.
    For wildcard segments over arrays, returns one pair per element.
    If any intermediate segment is missing, returns nothing (silent no-op).
    """
    if not segments:
        return []

    results = []

    def _walk(obj, depth):
        if depth == len(segments) - 1:
            seg = segments[depth]
            if seg == "*":
                if isinstance(obj, list):
                    for i in range(len(obj)):
                        results.append((obj, i))
            elif isinstance(obj, dict) and seg in obj:
                results.append((obj, seg))
            return

        seg = segments[depth]
        if seg == "*":
            if isinstance(obj, list):
                for item in obj:
                    _walk(item, depth + 1)
        elif isinstance(obj, dict) and seg in obj:
            _walk(obj[seg], depth + 1)

    _walk(data, 0)
    return results


def _ensure_parents(data, segments):
    """Walk the path, creating intermediate dicts as needed. Return (parent, final_key)."""
    obj = data
    for seg in segments[:-1]:
        if seg not in obj or not isinstance(obj[seg], dict):
            obj[seg] = {}
        obj = obj[seg]
    return obj, segments[-1]


def _apply_set(data, op):
    segments = _parse_path(op["path"])
    value = op["value"]
    if "*" not in segments:
        parent, key = _ensure_parents(data, segments)
        parent[key] = value
    else:
        _walk_and_set(data, segments, 0, value)


def _walk_and_set(obj, segments, depth, value):
    """Recursively walk segments, expanding wildcards over lists,
    and set the value at the final segment (creating it if absent)."""
    if depth == len(segments) - 1:
        seg = segments[depth]
        if seg == "*":
            if isinstance(obj, list):
                for i in range(len(obj)):
                    obj[i] = value
        elif isinstance(obj, dict):
            obj[seg] = value
        return

    seg = segments[depth]
    if seg == "*":
        if isinstance(obj, list):
            for item in obj:
                _walk_and_set(item, segments, depth + 1, value)
    elif isinstance(obj, dict):
        # With a wildcard elsewhere in the path, a missing intermediate is a
        # silent no-op. Don't materialize empty containers that would then
        # fail to match the wildcard anyway.
        if seg in obj and isinstance(obj[seg], (dict, list)):
            _walk_and_set(obj[seg], segments, depth + 1, value)


def _apply_remove(data, op):
    segments = _parse_path(op["path"])
    for parent, key in _resolve(data, segments):
        if isinstance(parent, dict):
            parent.pop(key, None)


def _apply_rename_key(data, op):
    segments = _parse_path(op["path"])
    to = op["to"]
    for parent, key in _resolve(data, segments):
        if isinstance(parent, dict) and key in parent:
            if to in parent:
                raise ValueError(
                    "rename_key: target key '" + to + "' already exists at path '" + op["path"] + "'"
                )
            parent[to] = parent.pop(key)


def _apply_map(data, op):
    segments = _parse_path(op["path"])
    from_val = op["from"]
    to_val = op["to"]
    for parent, key in _resolve(data, segments):
        if parent[key] == from_val:
            parent[key] = to_val


def _apply_move(data, op):
    segments = _parse_path(op["path"])
    matches = _resolve(data, segments)
    if not matches:
        return
    parent, key = matches[0]
    value = parent.pop(key)
    to_segments = _parse_path(op["to"])
    dest_parent, dest_key = _ensure_parents(data, to_segments)
    dest_parent[dest_key] = value


_OPS = {
    "set": _apply_set,
    "remove": _apply_remove,
    "rename_key": _apply_rename_key,
    "map": _apply_map,
    "move": _apply_move,
}


def _get_version(data):
    """Return the current schema_version, or None if absent."""
    metadata = data.get("metadata", {})
    if not isinstance(metadata, dict):
        return None
    return metadata.get("schema_version", None)


def _load_migrations(path=None):
    if path is None:
        path = _MIGRATIONS_PATH
    with open(path, "r") as f:
        return json.load(f)


def update_to_latest_schema(data, migrations_path=None):
    """Apply all applicable migration blocks to bring data to the latest schema version.

    Modifies data in place and returns it.
    """
    migrations = _load_migrations(migrations_path)

    while True:
        version = _get_version(data)
        applied = False
        for block in migrations:
            if block["from_version"] == version:
                for op in block["operations"]:
                    handler = _OPS.get(op["op"])
                    if handler is None:
                        raise ValueError("Unknown operation: " + op["op"])
                    handler(data, op)
                applied = True
                break
        if not applied:
            break

    return data