function tests = testMigrate
% Unit tests for migration-code/load_sample.m.
%
% These tests are intentionally kept outside GitHub Actions — MATLAB is
% proprietary and CI runners are not generally available. Run locally:
%
%   >> cd tests/matlab
%   >> results = runtests('testMigrate')
%
% Exercises the same invariants as the Python/Julia/JS suites so every
% language stays in lock-step.

    tests = functiontests(localfunctions);
end


% ── fixture helpers ───────────────────────────────────────────────────────

function repoRoot = getRepoRoot()
    thisDir = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(thisDir, '..', '..');
end

function setupOnce(testCase)
    repoRoot = getRepoRoot();
    addpath(fullfile(repoRoot, 'migration-code'));
    testCase.TestData.repoRoot = repoRoot;
    testCase.TestData.fixturesDir = fullfile(repoRoot, 'tests', 'fixtures');
    testCase.TestData.patchPath = fullfile(repoRoot, 'current', 'patch.json');
end

function data = loadFixture(testCase, name)
    data = load_sample(...
        fullfile(testCase.TestData.fixturesDir, name), ...
        testCase.TestData.patchPath);
end

function n = arrayLen(arr)
    if iscell(arr)
        n = numel(arr);
    elseif isstruct(arr)
        n = numel(arr);
    else
        n = 0;
    end
end

function elem = arrayGet(arr, idx)
    if iscell(arr)
        elem = arr{idx};
    else
        elem = arr(idx);
    end
end


% ── invariants ────────────────────────────────────────────────────────────

function verifyCurrent(testCase, data)
    verifyEqual(testCase, data.metadata.schema_version, '0.4.0');
    verifyEqual(testCase, data.metadata.schema_source, ...
        'https://raw.githubusercontent.com/nmr-samples/schema/main/versions/v0.4.0/schema.json');
end


% ── tests ─────────────────────────────────────────────────────────────────

function testV002MultiComponentRenamesEveryKey(testCase)
    data = loadFixture(testCase, 'sample_v0.0.2_multi.json');
    verifyCurrent(testCase, data);

    verifyFalse(testCase, isfield(data, 'Sample'));
    verifyFalse(testCase, isfield(data, 'Buffer'));
    verifyFalse(testCase, isfield(data, 'NMR_Tube'));
    verifyTrue(testCase, isfield(data, 'sample'));
    verifyTrue(testCase, isfield(data, 'people'));

    comps = data.sample.components;
    verifyEqual(testCase, arrayLen(comps), 3);
    for k = 1:arrayLen(comps)
        c = arrayGet(comps, k);
        verifyFalse(testCase, isfield(c, 'Name'));
        verifyFalse(testCase, isfield(c, 'Concentration'));
        verifyFalse(testCase, isfield(c, 'Unit'));
        verifyFalse(testCase, isfield(c, 'Isotopic_labelling'));
        verifyTrue(testCase, isfield(c, 'name'));
        verifyTrue(testCase, isfield(c, 'concentration_or_amount'));
        verifyTrue(testCase, isfield(c, 'unit'));
        verifyTrue(testCase, isfield(c, 'isotopic_labelling'));
        verifyTrue(testCase, isfield(c, 'molecular_weight'));
        verifyTrue(testCase, isfield(c, 'type'));
        verifyTrue(testCase, isempty(c.type));
        verifyNotEqual(testCase, c.unit, 'equiv');
    end

    verifyFalse(testCase, isfield(data.nmr_tube, 'diameter'));
    verifyEqual(testCase, data.nmr_tube.diameter_mm, 5.0);
end


function testV020WildcardMapAndSetHitEveryElement(testCase)
    data = loadFixture(testCase, 'sample_v0.2.0_multi.json');
    verifyCurrent(testCase, data);

    comps = data.sample.components;
    verifyEqual(testCase, arrayLen(comps), 3);
    for k = 1:arrayLen(comps)
        c = arrayGet(comps, k);
        verifyNotEqual(testCase, c.unit, 'equiv');
        verifyTrue(testCase, isfield(c, 'molecular_weight'));
        verifyTrue(testCase, isfield(c, 'type'));
        verifyTrue(testCase, isempty(c.type));
    end
end


function testV030AddsTypeToEveryComponent(testCase)
    data = loadFixture(testCase, 'sample_v0.3.0_multi.json');
    verifyCurrent(testCase, data);

    comps = data.sample.components;
    verifyEqual(testCase, arrayLen(comps), 3);
    for k = 1:arrayLen(comps)
        c = arrayGet(comps, k);
        verifyTrue(testCase, isfield(c, 'type'));
        verifyTrue(testCase, isempty(c.type));
        verifyTrue(testCase, isfield(c, 'molecular_weight'));
        verifyTrue(testCase, isfield(c, 'isotopic_labelling'));
    end
end


function testEmptyComponentsArrayIsNoopForWildcardOps(testCase)
    data = loadFixture(testCase, 'sample_v0.2.0_empty_components.json');
    verifyCurrent(testCase, data);
    verifyTrue(testCase, isempty(data.sample.components));
end


function testMissingComponentsIsSilentNoop(testCase)
    data = loadFixture(testCase, 'sample_v0.2.0_no_components.json');
    verifyCurrent(testCase, data);
    verifyFalse(testCase, isfield(data.sample, 'components'));
end


function testAlreadyCurrentIsNoop(testCase)
    data = loadFixture(testCase, 'sample_v0.4.0_already_current.json');
    verifyCurrent(testCase, data);
    % single-component JSON arrays decode as a struct in MATLAB
    if iscell(data.sample.components)
        c = data.sample.components{1};
    else
        c = data.sample.components(1);
    end
    verifyEqual(testCase, c.isotopic_labelling, '19F');
    verifyEqual(testCase, c.type, 'protein');
end


function testMigrationIsIdempotent(testCase)
    first  = loadFixture(testCase, 'sample_v0.0.2_multi.json');
    second = loadFixture(testCase, 'sample_v0.0.2_multi.json');
    verifyEqual(testCase, first.metadata.schema_version, '0.4.0');
    verifyEqual(testCase, arrayLen(first.sample.components), 3);
    verifyEqual(testCase, first.metadata.schema_version, second.metadata.schema_version);
end
