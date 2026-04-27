"""
Migrate NMR sample metadata to the latest schema version.

Call loadsample(path, migrations_path=MIGRATIONS_PATH[]) to load a JSON file
and apply migrations. It returns the migrated data as a Dict.

Call updatetolatestschema!(data) with a parsed JSON Dict. It modifies
the dict in place and returns it.

The migration patch file is expected at current/patch.json relative to
this module. Override by setting MIGRATIONS_PATH[] or passing migrations_path.
"""
module SchemaMigrate

using JSON

const MIGRATIONS_PATH = Ref(joinpath(@__DIR__, "..", "current", "patch.json"))


function _parse_path(path)
    isempty(path) && return String[]
    startswith(path, "/") || error("Path must start with '/': " * path)
    parts = split(path[2:end], "/")
    return [replace(replace(p, "~1" => "/"), "~0" => "~") for p in parts]
end


function _resolve(data, segments)
    results = Tuple{Any,Any}[]
    isempty(segments) && return results

    function _walk(obj, depth)
        if depth == length(segments)
            seg = segments[depth]
            if seg == "*"
                if isa(obj, Vector)
                    for i in eachindex(obj)
                        push!(results, (obj, i))
                    end
                end
            elseif isa(obj, Dict) && haskey(obj, seg)
                push!(results, (obj, seg))
            end
            return
        end

        seg = segments[depth]
        if seg == "*"
            if isa(obj, Vector)
                for item in obj
                    _walk(item, depth + 1)
                end
            end
        elseif isa(obj, Dict) && haskey(obj, seg)
            _walk(obj[seg], depth + 1)
        end
    end

    _walk(data, 1)
    return results
end


function _ensure_parents(data, segments)
    obj = data
    for seg in segments[1:end-1]
        if !haskey(obj, seg) || !isa(obj[seg], Dict)
            obj[seg] = Dict{String,Any}()
        end
        obj = obj[seg]
    end
    return obj, segments[end]
end


function _apply_set(data, op)
    segments = _parse_path(op["path"])
    value = op["value"]
    if !("*" in segments)
        parent, key = _ensure_parents(data, segments)
        parent[key] = value
    else
        _walk_and_set(data, segments, 1, value)
    end
end


function _walk_and_set(obj, segments, depth, value)
    if depth == length(segments)
        seg = segments[depth]
        if seg == "*"
            if isa(obj, Vector)
                for i in eachindex(obj)
                    obj[i] = value
                end
            end
        elseif isa(obj, Dict)
            obj[seg] = value
        end
        return
    end

    seg = segments[depth]
    if seg == "*"
        if isa(obj, Vector)
            for item in obj
                _walk_and_set(item, segments, depth + 1, value)
            end
        end
    elseif isa(obj, Dict)
        # With a wildcard elsewhere in the path, a missing intermediate is a
        # silent no-op. Don't materialize empty containers.
        if haskey(obj, seg) && isa(obj[seg], Union{Dict,Vector})
            _walk_and_set(obj[seg], segments, depth + 1, value)
        end
    end
end


function _apply_remove(data, op)
    segments = _parse_path(op["path"])
    for (parent, key) in _resolve(data, segments)
        if isa(parent, Dict)
            delete!(parent, key)
        end
    end
end


function _apply_rename_key(data, op)
    segments = _parse_path(op["path"])
    to = op["to"]
    for (parent, key) in _resolve(data, segments)
        if isa(parent, Dict) && haskey(parent, key)
            if haskey(parent, to)
                error("rename_key: target key '" * to * "' already exists at path '" * op["path"] * "'")
            end
            parent[to] = pop!(parent, key)
        end
    end
end


function _apply_map(data, op)
    segments = _parse_path(op["path"])
    from_val = op["from"]
    to_val = op["to"]
    for (parent, key) in _resolve(data, segments)
        if parent[key] == from_val
            parent[key] = to_val
        end
    end
end


function _apply_move(data, op)
    segments = _parse_path(op["path"])
    matches = _resolve(data, segments)
    isempty(matches) && return
    parent, key = matches[1]
    value = pop!(parent, key)
    to_segments = _parse_path(op["to"])
    dest_parent, dest_key = _ensure_parents(data, to_segments)
    dest_parent[dest_key] = value
end


const _OPS = Dict{String,Function}(
    "set" => _apply_set,
    "remove" => _apply_remove,
    "rename_key" => _apply_rename_key,
    "map" => _apply_map,
    "move" => _apply_move,
)


function _get_version(data)
    metadata = get(data, "metadata", nothing)
    isa(metadata, Dict) || return nothing
    return get(metadata, "schema_version", nothing)
end


function _load_migrations(path=MIGRATIONS_PATH[])
    return JSON.parsefile(path)
end

"""
    updatetolatestschema!(data, migrations_path=MIGRATIONS_PATH[]) -> Dict

Apply migrations to the given data Dict to update it to the latest schema version.
The data is modified in place and returned. Migrations are loaded from the given
migrations_path (default is MIGRATIONS_PATH[]).
"""
function updatetolatestschema!(data, migrations_path=MIGRATIONS_PATH[])
    migrations = _load_migrations(migrations_path)

    while true
        version = _get_version(data)
        applied = false
        for block in migrations
            if block["from_version"] == version
                for op in block["operations"]
                    handler = get(_OPS, op["op"], nothing)
                    handler === nothing && error("Unknown operation: " * op["op"])
                    handler(data, op)
                end
                applied = true
                break
            end
        end
        applied || break
    end

    return data
end

"""
    loadsample(path, migrations_path=MIGRATIONS_PATH[]) -> Dict

Load a JSON file from the given path and apply migrations to update it to the
latest schema version. Returns the migrated data as a Dict.
"""
function loadsample(path, migrations_path=MIGRATIONS_PATH[])
    data = JSON.parsefile(path; dicttype=Dict{String,Any})
    return updatetolatestschema!(data, migrations_path)
end

export updatetolatestschema!, loadsample

end # module