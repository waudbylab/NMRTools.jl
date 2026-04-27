/**
 * Migrate NMR sample metadata to the latest schema version.
 *
 * Call updateToLatestSchema(data, migrations) with a parsed JSON object
 * and a loaded migrations array. It modifies the object in place and
 * returns it. This function is synchronous.
 *
 * Call loadMigrations() to fetch the migrations array. By default this
 * fetches from GitHub; pass a URL to override.
 *
 * Call loadSample(migrations) to open a sample JSON file via the File
 * System Access API. Returns the parsed and migrated object.
 */

var MIGRATIONS_URL = "https://raw.githubusercontent.com/nmr-samples/schema/main/current/patch.json";


// ── Core migration engine (synchronous, no I/O) ─────────────────────

function _parsePath(path) {
    if (!path) return [];
    if (path[0] !== "/") throw new Error("Path must start with '/': " + path);
    var parts = path.slice(1).split("/");
    return parts.map(function (p) {
        return p.replace(/~1/g, "/").replace(/~0/g, "~");
    });
}


function _resolve(data, segments) {
    var results = [];
    if (!segments.length) return results;

    function _walk(obj, depth) {
        if (depth === segments.length - 1) {
            var seg = segments[depth];
            if (seg === "*") {
                if (Array.isArray(obj)) {
                    for (var i = 0; i < obj.length; i++) {
                        results.push([obj, i]);
                    }
                }
            } else if (obj !== null && typeof obj === "object" && !Array.isArray(obj) && seg in obj) {
                results.push([obj, seg]);
            }
            return;
        }

        var seg = segments[depth];
        if (seg === "*") {
            if (Array.isArray(obj)) {
                for (var i = 0; i < obj.length; i++) {
                    _walk(obj[i], depth + 1);
                }
            }
        } else if (obj !== null && typeof obj === "object" && !Array.isArray(obj) && seg in obj) {
            _walk(obj[seg], depth + 1);
        }
    }

    _walk(data, 0);
    return results;
}


function _ensureParents(data, segments) {
    var obj = data;
    for (var i = 0; i < segments.length - 1; i++) {
        var seg = segments[i];
        if (!(seg in obj) || typeof obj[seg] !== "object" || obj[seg] === null || Array.isArray(obj[seg])) {
            obj[seg] = {};
        }
        obj = obj[seg];
    }
    return [obj, segments[segments.length - 1]];
}


function _applySet(data, op) {
    var segments = _parsePath(op.path);
    var value = op.value;
    if (segments.indexOf("*") === -1) {
        var pair = _ensureParents(data, segments);
        pair[0][pair[1]] = value;
    } else {
        _walkAndSet(data, segments, 0, value);
    }
}


function _walkAndSet(obj, segments, depth, value) {
    if (depth === segments.length - 1) {
        var seg = segments[depth];
        if (seg === "*") {
            if (Array.isArray(obj)) {
                for (var i = 0; i < obj.length; i++) {
                    obj[i] = value;
                }
            }
        } else if (obj !== null && typeof obj === "object" && !Array.isArray(obj)) {
            obj[seg] = value;
        }
        return;
    }

    var seg = segments[depth];
    if (seg === "*") {
        if (Array.isArray(obj)) {
            for (var i = 0; i < obj.length; i++) {
                _walkAndSet(obj[i], segments, depth + 1, value);
            }
        }
    } else if (obj !== null && typeof obj === "object" && !Array.isArray(obj)) {
        // With a wildcard elsewhere in the path, a missing intermediate is
        // a silent no-op. Don't materialize empty containers.
        if (seg in obj && obj[seg] !== null && typeof obj[seg] === "object") {
            _walkAndSet(obj[seg], segments, depth + 1, value);
        }
    }
}


function _applyRemove(data, op) {
    var segments = _parsePath(op.path);
    var matches = _resolve(data, segments);
    for (var i = 0; i < matches.length; i++) {
        var parent = matches[i][0];
        var key = matches[i][1];
        if (typeof parent === "object" && !Array.isArray(parent)) {
            delete parent[key];
        }
    }
}


function _applyRenameKey(data, op) {
    var segments = _parsePath(op.path);
    var to = op.to;
    var matches = _resolve(data, segments);
    for (var i = 0; i < matches.length; i++) {
        var parent = matches[i][0];
        var key = matches[i][1];
        if (typeof parent === "object" && !Array.isArray(parent) && key in parent) {
            if (to in parent) {
                throw new Error("rename_key: target key '" + to + "' already exists at path '" + op.path + "'");
            }
            parent[to] = parent[key];
            delete parent[key];
        }
    }
}


function _applyMap(data, op) {
    var segments = _parsePath(op.path);
    var fromVal = op.from;
    var toVal = op.to;
    var matches = _resolve(data, segments);
    for (var i = 0; i < matches.length; i++) {
        var parent = matches[i][0];
        var key = matches[i][1];
        if (parent[key] === fromVal) {
            parent[key] = toVal;
        }
    }
}


function _applyMove(data, op) {
    var segments = _parsePath(op.path);
    var matches = _resolve(data, segments);
    if (!matches.length) return;
    var parent = matches[0][0];
    var key = matches[0][1];
    var value = parent[key];
    delete parent[key];
    var toSegments = _parsePath(op.to);
    var pair = _ensureParents(data, toSegments);
    pair[0][pair[1]] = value;
}


var _OPS = {
    "set": _applySet,
    "remove": _applyRemove,
    "rename_key": _applyRenameKey,
    "map": _applyMap,
    "move": _applyMove
};


function _getVersion(data) {
    if (!data.metadata || typeof data.metadata !== "object") return null;
    var v = data.metadata.schema_version;
    return v === undefined ? null : v;
}


/**
 * Apply all applicable migrations to data. Synchronous.
 * @param {Object} data - Parsed JSON sample object (modified in place)
 * @param {Array} migrations - Loaded migrations array (from loadMigrations)
 * @returns {Object} The migrated data
 */
function updateToLatestSchema(data, migrations) {
    while (true) {
        var version = _getVersion(data);
        var applied = false;
        for (var i = 0; i < migrations.length; i++) {
            var block = migrations[i];
            if (block.from_version === version) {
                for (var j = 0; j < block.operations.length; j++) {
                    var op = block.operations[j];
                    var handler = _OPS[op.op];
                    if (!handler) throw new Error("Unknown operation: " + op.op);
                    handler(data, op);
                }
                applied = true;
                break;
            }
        }
        if (!applied) break;
    }
    return data;
}


/**
 * Fetch the migrations array from a URL.
 * @param {string} [url] - Defaults to MIGRATIONS_URL (GitHub raw)
 * @returns {Promise<Array>}
 */
async function loadMigrations(url) {
    var response = await fetch(url || MIGRATIONS_URL);
    if (!response.ok) throw new Error("Failed to fetch migrations: " + response.status);
    return response.json();
}


/**
 * Open a sample JSON file via the File System Access API, parse and migrate it.
 * @param {Array} migrations - Loaded migrations array (from loadMigrations)
 * @returns {Promise<Object>} The migrated data
 */
async function loadSample(migrations) {
    var picks = await window.showOpenFilePicker({
        types: [{ description: "JSON", accept: { "application/json": [".json"] } }]
    });
    var file = await picks[0].getFile();
    var text = await file.text();
    var data = JSON.parse(text);
    return updateToLatestSchema(data, migrations);
}


// Node.js export (harmless in browsers: `module` is undefined there).
if (typeof module !== "undefined" && module.exports) {
    module.exports = {
        MIGRATIONS_URL: MIGRATIONS_URL,
        updateToLatestSchema: updateToLatestSchema,
        loadMigrations: loadMigrations,
        loadSample: loadSample,
        _parsePath: _parsePath,
        _resolve: _resolve,
        _applySet: _applySet,
        _applyRemove: _applyRemove,
        _applyRenameKey: _applyRenameKey,
        _applyMap: _applyMap,
        _applyMove: _applyMove
    };
}
