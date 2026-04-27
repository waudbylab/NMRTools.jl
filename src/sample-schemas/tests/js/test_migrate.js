"use strict";

const { test } = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

const REPO_ROOT = path.resolve(__dirname, "..", "..");
const FIXTURES_DIR = path.join(REPO_ROOT, "tests", "fixtures");
const PATCH_PATH = path.join(REPO_ROOT, "current", "patch.json");

const migrate = require(path.join(REPO_ROOT, "migration-code", "schema_migrate.js"));

const MIGRATIONS = JSON.parse(fs.readFileSync(PATCH_PATH, "utf8"));

const LATEST_VERSION = "0.4.0";
const LATEST_SCHEMA_SOURCE =
    "https://raw.githubusercontent.com/nmr-samples/schema/main/versions/v0.4.0/schema.json";

function load(name) {
    return JSON.parse(fs.readFileSync(path.join(FIXTURES_DIR, name), "utf8"));
}

function assertCurrent(data) {
    assert.equal(data.metadata.schema_version, LATEST_VERSION);
    assert.equal(data.metadata.schema_source, LATEST_SCHEMA_SOURCE);
}

function clone(obj) {
    return JSON.parse(JSON.stringify(obj));
}

// ── Core DSL unit tests ─────────────────────────────────────────────────────

test("set creates intermediate dicts on concrete paths", () => {
    const data = {};
    migrate._applySet(data, { op: "set", path: "/a/b/c", value: 42 });
    assert.deepEqual(data, { a: { b: { c: 42 } } });
});

test("set wildcard on array adds field to every element", () => {
    const data = { items: [{ name: "a" }, { name: "b" }, { name: "c" }] };
    migrate._applySet(data, { op: "set", path: "/items/*/flag", value: null });
    assert.deepEqual(data.items, [
        { name: "a", flag: null },
        { name: "b", flag: null },
        { name: "c", flag: null },
    ]);
});

test("set wildcard on empty array is no-op", () => {
    const data = { items: [] };
    migrate._applySet(data, { op: "set", path: "/items/*/flag", value: null });
    assert.deepEqual(data, { items: [] });
});

test("set wildcard on missing key is silent no-op", () => {
    const data = {};
    migrate._applySet(data, { op: "set", path: "/items/*/flag", value: null });
    assert.deepEqual(data, {});
});

test("rename_key wildcard renames on every element", () => {
    const data = { items: [{ Old: 1 }, { Old: 2 }, { Other: 3 }] };
    migrate._applyRenameKey(
        data,
        { op: "rename_key", path: "/items/*/Old", to: "new" }
    );
    assert.deepEqual(data.items, [{ new: 1 }, { new: 2 }, { Other: 3 }]);
});

test("map wildcard replaces only matching values", () => {
    const data = { items: [{ u: "equiv" }, { u: "mM" }, { u: "equiv" }] };
    migrate._applyMap(
        data,
        { op: "map", path: "/items/*/u", from: "equiv", to: "" }
    );
    assert.deepEqual(data.items, [{ u: "" }, { u: "mM" }, { u: "" }]);
});

test("remove wildcard drops key from each element", () => {
    const data = { items: [{ keep: 1, drop: 2 }, { keep: 3, drop: 4 }] };
    migrate._applyRemove(data, { op: "remove", path: "/items/*/drop" });
    assert.deepEqual(data.items, [{ keep: 1 }, { keep: 3 }]);
});

test("move relocates value", () => {
    const data = { old: { k: 1 } };
    migrate._applyMove(data, { op: "move", path: "/old", to: "/new/inner" });
    assert.deepEqual(data, { new: { inner: { k: 1 } } });
});

test("parse_path rejects missing leading slash", () => {
    assert.throws(() => migrate._parsePath("no-leading-slash"));
});

test("parse_path unescapes pointer syntax", () => {
    assert.deepEqual(migrate._parsePath("/a~1b/c~0d"), ["a/b", "c~d"]);
});

// ── End-to-end migration tests ──────────────────────────────────────────────

test("v0.0.2 multi-component: renames every key on every array entry", () => {
    const data = load("sample_v0.0.2_multi.json");
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assertCurrent(data);

    for (const key of ["Sample", "Buffer", "NMR Tube", "Laboratory Reference", "Notes"]) {
        assert.ok(!(key in data), `unexpected key ${key} after migration`);
    }
    assert.ok("sample" in data);
    assert.ok("people" in data);
    assert.deepEqual(data.people.users, ["Alice", "Bob"]);

    const comps = data.sample.components;
    assert.equal(comps.length, 3);
    for (const c of comps) {
        assert.ok(!("Name" in c));
        assert.ok(!("Concentration" in c));
        assert.ok(!("Unit" in c));
        assert.ok(!("Isotopic labelling" in c));
        assert.ok(!("Custom labelling" in c));
        assert.ok("name" in c);
        assert.ok("concentration_or_amount" in c);
        assert.ok("unit" in c);
        assert.ok("isotopic_labelling" in c);
        assert.ok("custom_labelling" in c);
        assert.ok("molecular_weight" in c);
        assert.ok("type" in c);
        assert.equal(c.type, null);
    }
    for (const c of comps) assert.notEqual(c.unit, "equiv");

    // sample-level renames
    assert.ok(!("Label" in data.sample));
    assert.equal(data.sample.label, "Test v0.0.2 with multiple components");
    assert.equal(data.sample.physical_form, "");  // added in 0.0.3 → 0.1.0

    // buffer field renames
    assert.ok(!("pH" in data.buffer));
    assert.equal(data.buffer.ph, 7.4);
    assert.ok(!("Solvent" in data.buffer));
    assert.equal(data.buffer.solvent, "10% D2O");
    assert.ok(!("Chemical shift reference" in data.buffer));
    assert.equal(data.buffer.chemical_shift_reference, "DSS");
    assert.ok(!("Reference concentration" in data.buffer));
    assert.equal(data.buffer.reference_concentration, 10);
    assert.ok(!("Reference unit" in data.buffer));
    assert.equal(data.buffer.reference_unit, "uM");
    assert.ok(!("Custom solvent" in data.buffer));
    assert.equal(data.buffer.custom_solvent, "");

    // buffer component renames
    const bufComps = data.buffer.components;
    assert.equal(bufComps.length, 2);
    for (const bc of bufComps) {
        assert.ok(!("Concentration" in bc));
        assert.ok(!("Unit" in bc));
        assert.ok("concentration" in bc);
        assert.ok("unit" in bc);
    }

    // NMR tube renames
    assert.ok(!("Type" in data.nmr_tube));
    assert.equal(data.nmr_tube.type, "shigemi");
    assert.ok(!("Sample Volume (μL)" in data.nmr_tube));
    assert.equal(data.nmr_tube.sample_volume_uL, 300);
    assert.ok(!("samplejet_rack_position" in data.nmr_tube));  // removed in 0.0.3 → 0.1.0
    assert.ok(!("samplejet_rack_id" in data.nmr_tube));
    assert.equal(data.nmr_tube.rack_id, "rack-001");

    // diameter string mapped to number AND renamed to diameter_mm
    assert.ok(!("diameter" in data.nmr_tube));
    assert.equal(data.nmr_tube.diameter_mm, 5.0);

    // reference field renames
    assert.ok(!("Labbook Entry" in data.reference));
    assert.equal(data.reference.labbook_entry, "page 42");
    assert.ok(!("Experiment ID" in data.reference));
    assert.equal(data.reference.sample_id, "EXP-2024-001");
});

test("v0.2.0 multi-component: wildcard map+set hit every component", () => {
    const data = load("sample_v0.2.0_multi.json");
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assertCurrent(data);

    const comps = data.sample.components;
    assert.equal(comps.length, 3);
    for (const c of comps) {
        assert.notEqual(c.unit, "equiv");
        assert.ok("molecular_weight" in c);
        assert.ok("type" in c);
        assert.equal(c.type, null);
    }

    // nmr_tube.diameter renamed to diameter_mm in 0.2.0 → 0.3.0
    assert.ok(!("diameter" in data.nmr_tube));
    assert.ok("diameter_mm" in data.nmr_tube);
    assert.equal(data.nmr_tube.diameter_mm, 5.0);
});

test("v0.3.0 multi-component: adds type to every component", () => {
    const data = load("sample_v0.3.0_multi.json");
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assertCurrent(data);

    const comps = data.sample.components;
    assert.equal(comps.length, 3);
    for (const c of comps) {
        assert.ok("type" in c);
        assert.equal(c.type, null);
        assert.ok("molecular_weight" in c);
    }
});

for (const [diameterStr, expected] of [["1.7 mm", 1.7], ["3 mm", 3.0], ["5 mm", 5.0], ["", null]]) {
    test(`v0.0.3 diameter "${diameterStr}" maps to ${expected}`, () => {
        const data = {
            sample: {},
            nmr_tube: { diameter: diameterStr },
            metadata: { schema_version: "0.0.3" },
        };
        migrate.updateToLatestSchema(data, MIGRATIONS);
        assert.ok(!("diameter" in data.nmr_tube));
        assert.equal(data.nmr_tube.diameter_mm, expected);
    });
}

test("v0.0.3 samplejet rename, remove, and physical_form addition", () => {
    const data = {
        sample: {},
        nmr_tube: {
            diameter: "5 mm",
            samplejet_rack_id: "rack-001",
            samplejet_rack_position: "A3",
        },
        metadata: { schema_version: "0.0.3" },
    };
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assert.equal(data.nmr_tube.rack_id, "rack-001");
    assert.ok(!("samplejet_rack_id" in data.nmr_tube));
    assert.ok(!("samplejet_rack_position" in data.nmr_tube));
    assert.equal(data.sample.physical_form, "");
});

test("v0.1.0 chain migrates to current", () => {
    const data = {
        sample: {
            physical_form: "solution",
            components: [{ concentration_or_amount: 1.0, unit: "mM", isotopic_labelling: "13C,15N" }],
        },
        nmr_tube: { diameter: 5.0 },
        metadata: { schema_version: "0.1.0" },
    };
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assertCurrent(data);
    assert.ok(!("diameter" in data.nmr_tube));
    assert.equal(data.nmr_tube.diameter_mm, 5.0);
});

test("empty components array: wildcard ops are no-ops", () => {
    const data = load("sample_v0.2.0_empty_components.json");
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assertCurrent(data);
    assert.deepEqual(data.sample.components, []);
});

test("missing components key: wildcard ops are silent no-ops", () => {
    const data = load("sample_v0.2.0_no_components.json");
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assertCurrent(data);
    assert.ok(!("components" in data.sample));
});

test("already-current sample is unchanged", () => {
    const data = load("sample_v0.4.0_already_current.json");
    const before = clone(data);
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assert.deepEqual(data, before);
    assert.equal(data.sample.components[0].isotopic_labelling, "19F");
});

test("migration is idempotent", () => {
    const data = load("sample_v0.0.2_multi.json");
    migrate.updateToLatestSchema(data, MIGRATIONS);
    const first = clone(data);
    migrate.updateToLatestSchema(data, MIGRATIONS);
    assert.deepEqual(data, first);
});

test("stress: 50-element component array", () => {
    const data = {
        sample: {
            components: Array.from({ length: 50 }, (_, i) => ({
                Name: `c${i}`,
                Concentration: i,
                Unit: i % 2 === 0 ? "mM" : "equiv",
            })),
        },
        metadata: { schema_version: "0.0.2" },
    };
    migrate.updateToLatestSchema(data, MIGRATIONS);
    const comps = data.sample.components;
    assert.equal(comps.length, 50);
    for (const c of comps) {
        assert.ok("name" in c);
        assert.ok("concentration_or_amount" in c);
        assert.notEqual(c.unit, "equiv");
        assert.equal(c.type, null);
        assert.ok("molecular_weight" in c);
    }
});
