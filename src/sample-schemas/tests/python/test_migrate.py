"""Unit tests for the Python migration engine (migration-code/migrate.py)."""

import copy
import json
import os

import pytest

from conftest import FIXTURES_DIR, PATCH_PATH

import migrate


LATEST_VERSION = "0.4.0"
LATEST_SCHEMA_SOURCE = (
    "https://raw.githubusercontent.com/nmr-samples/schema/main/versions/v0.4.0/schema.json"
)


def _load(name):
    with open(os.path.join(FIXTURES_DIR, name), "r") as f:
        return json.load(f)


def _migrate(data):
    return migrate.update_to_latest_schema(data, migrations_path=PATCH_PATH)


# ── Core DSL unit tests ─────────────────────────────────────────────────────


def test_set_creates_intermediate_dicts():
    data = {}
    migrate._apply_set(data, {"op": "set", "path": "/a/b/c", "value": 42})
    assert data == {"a": {"b": {"c": 42}}}


def test_set_wildcard_on_array_adds_field_to_every_element():
    data = {"items": [{"name": "a"}, {"name": "b"}, {"name": "c"}]}
    migrate._apply_set(data, {"op": "set", "path": "/items/*/flag", "value": None})
    assert data["items"] == [
        {"name": "a", "flag": None},
        {"name": "b", "flag": None},
        {"name": "c", "flag": None},
    ]


def test_set_wildcard_on_empty_array_is_noop():
    data = {"items": []}
    migrate._apply_set(data, {"op": "set", "path": "/items/*/flag", "value": None})
    assert data == {"items": []}


def test_set_wildcard_on_missing_key_is_noop():
    data = {}
    migrate._apply_set(data, {"op": "set", "path": "/items/*/flag", "value": None})
    assert data == {}


def test_rename_key_wildcard_renames_on_every_element():
    data = {"items": [{"Old": 1}, {"Old": 2}, {"Other": 3}]}
    migrate._apply_rename_key(
        data, {"op": "rename_key", "path": "/items/*/Old", "to": "new"}
    )
    assert data["items"] == [{"new": 1}, {"new": 2}, {"Other": 3}]


def test_map_wildcard_replaces_only_matching_values():
    data = {"items": [{"u": "equiv"}, {"u": "mM"}, {"u": "equiv"}]}
    migrate._apply_map(
        data, {"op": "map", "path": "/items/*/u", "from": "equiv", "to": ""}
    )
    assert data["items"] == [{"u": ""}, {"u": "mM"}, {"u": ""}]


def test_remove_wildcard_drops_key_from_each_element():
    data = {"items": [{"keep": 1, "drop": 2}, {"keep": 3, "drop": 4}]}
    migrate._apply_remove(data, {"op": "remove", "path": "/items/*/drop"})
    assert data["items"] == [{"keep": 1}, {"keep": 3}]


def test_move_relocates_value():
    data = {"old": {"k": 1}}
    migrate._apply_move(data, {"op": "move", "path": "/old", "to": "/new/inner"})
    assert data == {"new": {"inner": {"k": 1}}}


def test_parse_path_rejects_missing_leading_slash():
    with pytest.raises(ValueError):
        migrate._parse_path("no-leading-slash")


def test_parse_path_unescapes_pointer_syntax():
    assert migrate._parse_path("/a~1b/c~0d") == ["a/b", "c~d"]


# ── End-to-end migration tests ──────────────────────────────────────────────


def _assert_current(data):
    assert data["metadata"]["schema_version"] == LATEST_VERSION
    assert data["metadata"]["schema_source"] == LATEST_SCHEMA_SOURCE


def test_migrate_v002_multi_component_renames_every_key():
    data = _load("sample_v0.0.2_multi.json")
    _migrate(data)
    _assert_current(data)

    # top-level renames
    for old in ("Sample", "Buffer", "NMR Tube", "Laboratory Reference", "Notes", "Metadata", "Users"):
        assert old not in data
    assert "sample" in data and "buffer" in data and "nmr_tube" in data
    assert "people" in data and data["people"]["users"] == ["Alice", "Bob"]

    # per-component renames happen on ALL array entries
    comps = data["sample"]["components"]
    assert len(comps) == 3
    for c in comps:
        assert "Name" not in c and "Concentration" not in c and "Unit" not in c
        assert "Isotopic labelling" not in c
        assert "Custom labelling" not in c
        assert "name" in c
        assert "concentration_or_amount" in c
        assert "unit" in c
        assert "isotopic_labelling" in c
        assert "custom_labelling" in c
        # v0.3.0 added molecular_weight; v0.4.0 added type
        assert "molecular_weight" in c
        assert "type" in c and c["type"] is None

    # equiv unit is stripped across every component
    units = [c["unit"] for c in comps]
    assert "equiv" not in units

    # sample-level renames
    assert "Label" not in data["sample"]
    assert data["sample"]["label"] == "Test v0.0.2 with multiple components"
    assert data["sample"]["physical_form"] == ""  # added in 0.0.3 → 0.1.0

    # buffer field renames
    assert "pH" not in data["buffer"]
    assert data["buffer"]["ph"] == 7.4
    assert "Solvent" not in data["buffer"]
    assert data["buffer"]["solvent"] == "10% D2O"
    assert "Chemical shift reference" not in data["buffer"]
    assert data["buffer"]["chemical_shift_reference"] == "DSS"
    assert "Reference concentration" not in data["buffer"]
    assert data["buffer"]["reference_concentration"] == 10
    assert "Reference unit" not in data["buffer"]
    assert data["buffer"]["reference_unit"] == "uM"
    assert "Custom solvent" not in data["buffer"]
    assert data["buffer"]["custom_solvent"] == ""

    # buffer component renames (Concentration/Unit → concentration/unit)
    buf_comps = data["buffer"]["components"]
    assert len(buf_comps) == 2
    for bc in buf_comps:
        assert "Concentration" not in bc and "Unit" not in bc
        assert "concentration" in bc and "unit" in bc

    # NMR tube renames
    assert "Type" not in data["nmr_tube"]
    assert data["nmr_tube"]["type"] == "shigemi"
    assert "Sample Volume (μL)" not in data["nmr_tube"]
    assert data["nmr_tube"]["sample_volume_uL"] == 300
    assert "samplejet_rack_position" not in data["nmr_tube"]  # removed in 0.0.3 → 0.1.0
    assert "samplejet_rack_id" not in data["nmr_tube"]
    assert data["nmr_tube"]["rack_id"] == "rack-001"

    # diameter string gets mapped to number AND renamed to diameter_mm
    assert "diameter" not in data["nmr_tube"]
    assert data["nmr_tube"]["diameter_mm"] == 5.0

    # reference field renames
    assert "Labbook Entry" not in data["reference"]
    assert data["reference"]["labbook_entry"] == "page 42"
    assert "Experiment ID" not in data["reference"]
    assert data["reference"]["sample_id"] == "EXP-2024-001"


def test_migrate_v020_wildcard_map_and_set_hit_every_element():
    data = _load("sample_v0.2.0_multi.json")
    comps_before = copy.deepcopy(data["sample"]["components"])
    assert sum(1 for c in comps_before if c["unit"] == "equiv") == 2

    _migrate(data)
    _assert_current(data)

    comps = data["sample"]["components"]
    assert len(comps) == 3
    for c in comps:
        assert c["unit"] != "equiv"
        assert "molecular_weight" in c
        assert "type" in c and c["type"] is None

    # nmr_tube.diameter renamed to diameter_mm in 0.2.0 → 0.3.0
    assert "diameter" not in data["nmr_tube"]
    assert "diameter_mm" in data["nmr_tube"]
    assert data["nmr_tube"]["diameter_mm"] == 5.0


def test_migrate_v030_adds_type_to_every_component():
    data = _load("sample_v0.3.0_multi.json")
    _migrate(data)
    _assert_current(data)

    comps = data["sample"]["components"]
    assert len(comps) == 3
    for c in comps:
        assert "type" in c and c["type"] is None
        assert "molecular_weight" in c
        assert "isotopic_labelling" in c

    # "unlabelled" in the fixture must be remapped
    assert comps[1]["isotopic_labelling"] == "natural abundance"


def test_migrate_empty_components_is_noop_for_wildcard_ops():
    data = _load("sample_v0.2.0_empty_components.json")
    _migrate(data)
    _assert_current(data)
    assert data["sample"]["components"] == []


def test_migrate_missing_components_is_silent_noop():
    data = _load("sample_v0.2.0_no_components.json")
    _migrate(data)
    _assert_current(data)
    # no components key should have been invented by wildcard set
    assert "components" not in data["sample"]


def test_migrate_already_current_is_noop():
    data = _load("sample_v0.4.0_already_current.json")
    before = copy.deepcopy(data)
    _migrate(data)
    assert data == before


def test_migrate_is_idempotent():
    data = _load("sample_v0.0.2_multi.json")
    _migrate(data)
    first = copy.deepcopy(data)
    _migrate(data)
    assert data == first


def test_migrate_preserves_19f_labelling():
    data = _load("sample_v0.4.0_already_current.json")
    _migrate(data)
    assert data["sample"]["components"][0]["isotopic_labelling"] == "19F"


def test_migrate_v030_isotopic_labelling_mappings():
    data = _load("sample_v0.3.0_labelling_solvent.json")
    _migrate(data)
    _assert_current(data)

    comps = data["sample"]["components"]
    assert comps[0]["isotopic_labelling"] == "natural abundance"   # unlabelled → natural abundance
    assert comps[1]["isotopic_labelling"] == "Ile-13CH3,15N"       # Ile-δ1-13CH3,15N → Ile-13CH3,15N
    assert comps[2]["isotopic_labelling"] == "2H,Ile-13CH3"        # 2H,Ile-δ1-13CH3 → 2H,Ile-13CH3
    assert comps[3]["isotopic_labelling"] == ""                    # 2H,Leu/Val-13CH3 → "" (no equivalent)


def test_migrate_v030_solvent_d6_dmso():
    data = _load("sample_v0.3.0_labelling_solvent.json")
    _migrate(data)
    assert data["buffer"]["solvent"] == "DMSO-d6"


@pytest.mark.parametrize("diameter_str,expected", [
    ("1.7 mm", 1.7),
    ("3 mm", 3.0),
    ("5 mm", 5.0),
    ("", None),
])
def test_migrate_v003_diameter_string_to_number(diameter_str, expected):
    data = {
        "sample": {},
        "nmr_tube": {"diameter": diameter_str},
        "metadata": {"schema_version": "0.0.3"},
    }
    _migrate(data)
    assert "diameter" not in data["nmr_tube"]
    assert data["nmr_tube"]["diameter_mm"] == expected


def test_migrate_v003_samplejet_rename_remove_and_physical_form():
    data = {
        "sample": {},
        "nmr_tube": {
            "diameter": "5 mm",
            "samplejet_rack_id": "rack-001",
            "samplejet_rack_position": "A3",
        },
        "metadata": {"schema_version": "0.0.3"},
    }
    _migrate(data)
    assert data["nmr_tube"]["rack_id"] == "rack-001"
    assert "samplejet_rack_id" not in data["nmr_tube"]
    assert "samplejet_rack_position" not in data["nmr_tube"]
    assert data["sample"]["physical_form"] == ""


def test_migrate_v010_chain_migrates_to_current():
    """A v0.1.0 document traverses 0.1.0→0.2.0→0.3.0→0.4.0 cleanly."""
    data = {
        "sample": {
            "physical_form": "solution",
            "components": [
                {
                    "concentration_or_amount": 1.0,
                    "unit": "mM",
                    "isotopic_labelling": "13C,15N",
                },
            ],
        },
        "nmr_tube": {"diameter": 5.0},
        "metadata": {"schema_version": "0.1.0"},
    }
    _migrate(data)
    _assert_current(data)
    assert "diameter" not in data["nmr_tube"]
    assert data["nmr_tube"]["diameter_mm"] == 5.0


def test_migrate_v030_solvent_d4_methanol():
    data = {
        "sample": {"components": []},
        "buffer": {"solvent": "D4-methanol"},
        "metadata": {"schema_version": "0.3.0"},
    }
    _migrate(data)
    assert data["buffer"]["solvent"] == "Methanol-d4"


def test_migrate_v030_unaffected_labelling_values_preserved():
    """Values already valid in v0.4.0 must not be altered."""
    data = {
        "sample": {
            "components": [
                {"isotopic_labelling": "13C,15N"},
                {"isotopic_labelling": "2H,13C,15N"},
                {"isotopic_labelling": ""},
                {"isotopic_labelling": "custom"},
            ]
        },
        "metadata": {"schema_version": "0.3.0"},
    }
    _migrate(data)
    labels = [c["isotopic_labelling"] for c in data["sample"]["components"]]
    assert labels == ["13C,15N", "2H,13C,15N", "", "custom"]


def test_migrate_large_component_array():
    """Stress test: wildcard operations on a 50-element array."""
    data = {
        "sample": {
            "components": [
                {"Name": "c%d" % i, "Concentration": i, "Unit": "equiv" if i % 2 else "mM"}
                for i in range(50)
            ]
        },
        "metadata": {"schema_version": "0.0.2"},
    }
    _migrate(data)
    comps = data["sample"]["components"]
    assert len(comps) == 50
    for c in comps:
        assert "name" in c
        assert "concentration_or_amount" in c
        assert c["unit"] != "equiv"
        assert c["type"] is None
        assert "molecular_weight" in c
