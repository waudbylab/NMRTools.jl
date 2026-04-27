% NMR Sample Loading and Schema Migration - MATLAB Implementation
%
% Provides a function to load NMR sample JSON files and migrate them to the
% latest schema version using the patch.json migration rules.
%
% Usage:
%   data = load_sample('path/to/sample.json');
%   data = load_sample('path/to/sample.json', 'path/to/patch.json');

function data = load_sample(samplePath, patchPath)
    if nargin < 2
        thisDir = fileparts(mfilename('fullpath'));
        patchPath = fullfile(thisDir, '..', 'current', 'patch.json');
    end
    text = fileread(samplePath);
    data = jsondecode(text);
    migrations = jsondecode(fileread(patchPath));
    data = updateToLatestSchema(data, migrations);
end


function data = updateToLatestSchema(data, migrations)
    while true
        version = getVersion(data);
        applied = false;
        for i = 1:numel(migrations)
            block = getElement(migrations, i);
            if strcmp(block.from_version, version)
                ops = block.operations;
                for j = 1:numel(ops)
                    data = applyOp(data, getElement(ops, j));
                end
                applied = true;
                break;
            end
        end
        if ~applied
            break;
        end
    end
end


function elem = getElement(arr, idx)
    if iscell(arr)
        elem = arr{idx};
    else
        elem = arr(idx);
    end
end


% -------------------------------------------------------------------------
% Internal helpers
% -------------------------------------------------------------------------

function version = getVersion(data)
    version = '';
    if isstruct(data) && isfield(data, 'metadata')
        m = data.metadata;
        if isstruct(m) && isfield(m, 'schema_version')
            version = m.schema_version;
        end
    end
end


function segments = parsePath(path)
    if isempty(path)
        segments = {};
        return;
    end
    assert(path(1) == '/', 'Path must start with /: %s', path);
    parts = strsplit(path(2:end), '/');
    segments = strrep(strrep(parts, '~1', '/'), '~0', '~');
end


function data = applyOp(data, op)
    switch op.op
        case 'set'
            data = applySet(data, op);
        case 'remove'
            data = applyRemove(data, op);
        case 'rename_key'
            data = applyRenameKey(data, op);
        case 'map'
            data = applyMap(data, op);
        case 'move'
            data = applyMove(data, op);
        otherwise
            error('Unknown operation: %s', op.op);
    end
end


% --- set -----------------------------------------------------------------

function data = applySet(data, op)
    segs = parsePath(op.path);
    data = setAtPath(data, segs, op.value);
end

function obj = setAtPath(obj, segs, value)
    if isempty(segs)
        return;
    end
    key = segs{1};
    rest = segs(2:end);
    if strcmp(key, '*')
        if isempty(rest)
            % wildcard at leaf: replace each element
            if iscell(obj)
                for k = 1:numel(obj)
                    obj{k} = value;
                end
            end
        else
            % wildcard at intermediate level: recurse into each element
            if iscell(obj)
                for k = 1:numel(obj)
                    obj{k} = setAtPath(obj{k}, rest, value);
                end
            elseif isstruct(obj) && numel(obj) > 1
                cellArr = cell(1, numel(obj));
                for k = 1:numel(obj)
                    cellArr{k} = setAtPath(obj(k), rest, value);
                end
                obj = cellArr;
            end
        end
    elseif isempty(rest)
        obj = setField(obj, key, value);
    else
        child = getField(obj, key);
        if ~isstruct(child) && ~iscell(child)
            % With a wildcard elsewhere in the path, a missing intermediate
            % is a silent no-op. Otherwise create an empty struct so we can
            % descend into concrete paths like /metadata/schema_version.
            if any(strcmp(rest, '*'))
                return;
            end
            child = struct();
        end
        child = setAtPath(child, rest, value);
        obj = setField(obj, key, child);
    end
end


% --- remove --------------------------------------------------------------

function data = applyRemove(data, op)
    segs = parsePath(op.path);
    data = removeAtPath(data, segs);
end

function obj = removeAtPath(obj, segs)
    if isempty(segs) || ~isstruct(obj)
        return;
    end
    key = segs{1};
    safeKey = matlabKey(key);
    rest = segs(2:end);
    if strcmp(key, '*')
        return;
    end
    if ~isfield(obj, safeKey)
        return;
    end
    if isempty(rest)
        obj = rmfield(obj, safeKey);
    else
        child = obj.(safeKey);
        if strcmp(rest{1}, '*') && (iscell(child) || isstruct(child))
            child = applyToArray(child, rest(2:end), @removeAtPath);
        else
            child = removeAtPath(child, rest);
        end
        obj.(safeKey) = child;
    end
end


% --- rename_key ----------------------------------------------------------

function data = applyRenameKey(data, op)
    segs = parsePath(op.path);
    data = renameAtPath(data, segs, op.to);
end

function obj = renameAtPath(obj, segs, toKey)
    if isempty(segs) || ~isstruct(obj)
        return;
    end
    key = segs{1};
    safeKey = matlabKey(key);
    rest = segs(2:end);
    if strcmp(key, '*')
        return;
    end
    if ~isfield(obj, safeKey)
        return;
    end
    if isempty(rest)
        safeToKey = matlabKey(toKey);
        if isfield(obj, safeToKey)
            error('rename_key: target key ''%s'' already exists', toKey);
        end
        obj.(safeToKey) = obj.(safeKey);
        obj = rmfield(obj, safeKey);
    else
        child = obj.(safeKey);
        if strcmp(rest{1}, '*') && (iscell(child) || isstruct(child))
            child = applyToArray(child, rest(2:end), @(c,s) renameAtPath(c, s, toKey));
        else
            child = renameAtPath(child, rest, toKey);
        end
        obj.(safeKey) = child;
    end
end


% --- map -----------------------------------------------------------------

function data = applyMap(data, op)
    segs = parsePath(op.path);
    data = mapAtPath(data, segs, op.from, op.to);
end

function obj = mapAtPath(obj, segs, fromVal, toVal)
    if isempty(segs) || ~isstruct(obj)
        return;
    end
    key = segs{1};
    safeKey = matlabKey(key);
    rest = segs(2:end);
    if strcmp(key, '*')
        return;
    end
    if ~isfield(obj, safeKey)
        return;
    end
    if isempty(rest)
        current = obj.(safeKey);
        if isequal(current, fromVal)
            obj.(safeKey) = toVal;
        end
    else
        child = obj.(safeKey);
        if strcmp(rest{1}, '*') && (iscell(child) || isstruct(child))
            child = applyToArray(child, rest(2:end), @(c,s) mapAtPath(c, s, fromVal, toVal));
        else
            child = mapAtPath(child, rest, fromVal, toVal);
        end
        obj.(safeKey) = child;
    end
end


% --- move ----------------------------------------------------------------

function data = applyMove(data, op)
    fromSegs = parsePath(op.path);
    toSegs = parsePath(op.to);
    [data, value, found] = extractAtPath(data, fromSegs);
    if found
        data = setAtPath(data, toSegs, value);
    end
end

function [obj, value, found] = extractAtPath(obj, segs)
    found = false;
    value = [];
    if isempty(segs) || ~isstruct(obj)
        return;
    end
    key = segs{1};
    safeKey = matlabKey(key);
    rest = segs(2:end);
    if strcmp(key, '*') || ~isfield(obj, safeKey)
        return;
    end
    if isempty(rest)
        value = obj.(safeKey);
        found = true;
        obj = rmfield(obj, safeKey);
    else
        child = obj.(safeKey);
        [child, value, found] = extractAtPath(child, rest);
        obj.(safeKey) = child;
    end
end


% -------------------------------------------------------------------------
% Array helpers
% -------------------------------------------------------------------------

function arr = applyToArray(arr, segs, fn)
    if iscell(arr)
        for k = 1:numel(arr)
            arr{k} = fn(arr{k}, segs);
        end
    elseif isstruct(arr) && numel(arr) > 1
        cellArr = cell(1, numel(arr));
        for k = 1:numel(arr)
            cellArr{k} = fn(arr(k), segs);
        end
        arr = cellArr;
    else
        arr = fn(arr, segs);
    end
end


% -------------------------------------------------------------------------
% Field access helpers that handle jsondecode naming quirks
% -------------------------------------------------------------------------

function mk = matlabKey(key)
    % MATLAB jsondecode removes spaces using camelCase (capitalise char after
    % each space). Other invalid characters are replaced with underscores.
    mk = regexprep(key, ' ([a-zA-Z])', '${upper($1)}');
    mk = strrep(mk, ' ', '');
    mk = regexprep(mk, '[^a-zA-Z0-9_]', '_');
    if ~isempty(mk) && mk(1) >= '0' && mk(1) <= '9'
        mk = ['x' mk];
    end
end

function value = getField(obj, key)
    mk = matlabKey(key);
    if isstruct(obj) && isfield(obj, mk)
        value = obj.(mk);
    else
        value = [];
    end
end

function obj = setField(obj, key, value)
    mk = matlabKey(key);
    if ~isstruct(obj)
        obj = struct();
    end
    obj.(mk) = value;
end
