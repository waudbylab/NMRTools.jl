using Test
using JSON

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MIGRATION_CODE_DIR = joinpath(REPO_ROOT, "migration-code")
const FIXTURES_DIR = joinpath(REPO_ROOT, "tests", "fixtures")
const PATCH_PATH = joinpath(REPO_ROOT, "current", "patch.json")

include(joinpath(MIGRATION_CODE_DIR, "migrate.jl"))
using .SchemaMigrate

const LATEST_VERSION = "0.4.0"
const LATEST_SCHEMA_SOURCE = "https://raw.githubusercontent.com/nmr-samples/schema/main/versions/v0.4.0/schema.json"

load_fixture(name) = JSON.parsefile(joinpath(FIXTURES_DIR, name); dicttype=Dict{String,Any})
migrate!(data) = updatetolatestschema!(data, PATCH_PATH)

function assert_current(data)
    @test data["metadata"]["schema_version"] == LATEST_VERSION
    @test data["metadata"]["schema_source"] == LATEST_SCHEMA_SOURCE
end

@testset "SchemaMigrate" begin

    @testset "Core DSL — wildcards over arrays" begin
        # set wildcard on array adds field to every element
        data = Dict{String,Any}(
            "items" => Any[
                Dict{String,Any}("name" => "a"),
                Dict{String,Any}("name" => "b"),
                Dict{String,Any}("name" => "c"),
            ],
        )
        SchemaMigrate._apply_set(data, Dict("op" => "set", "path" => "/items/*/flag", "value" => nothing))
        @test all(haskey(c, "flag") && c["flag"] === nothing for c in data["items"])

        # set wildcard on empty array is no-op
        data = Dict{String,Any}("items" => Any[])
        SchemaMigrate._apply_set(data, Dict("op" => "set", "path" => "/items/*/flag", "value" => nothing))
        @test data["items"] == Any[]

        # set wildcard on missing key is silent no-op (no intermediate created)
        data = Dict{String,Any}()
        SchemaMigrate._apply_set(data, Dict("op" => "set", "path" => "/items/*/flag", "value" => nothing))
        @test !haskey(data, "items")

        # rename_key wildcard renames on every element
        data = Dict{String,Any}(
            "items" => Any[
                Dict{String,Any}("Old" => 1),
                Dict{String,Any}("Old" => 2),
                Dict{String,Any}("Other" => 3),
            ],
        )
        SchemaMigrate._apply_rename_key(
            data, Dict("op" => "rename_key", "path" => "/items/*/Old", "to" => "new"),
        )
        @test data["items"][1] == Dict{String,Any}("new" => 1)
        @test data["items"][2] == Dict{String,Any}("new" => 2)
        @test data["items"][3] == Dict{String,Any}("Other" => 3)

        # map wildcard replaces only matching values
        data = Dict{String,Any}(
            "items" => Any[
                Dict{String,Any}("u" => "equiv"),
                Dict{String,Any}("u" => "mM"),
                Dict{String,Any}("u" => "equiv"),
            ],
        )
        SchemaMigrate._apply_map(
            data, Dict("op" => "map", "path" => "/items/*/u", "from" => "equiv", "to" => ""),
        )
        @test [c["u"] for c in data["items"]] == ["", "mM", ""]
    end

    @testset "parse_path" begin
        @test_throws Exception SchemaMigrate._parse_path("no-leading-slash")
        @test SchemaMigrate._parse_path("/a~1b/c~0d") == ["a/b", "c~d"]
    end

    @testset "v0.0.2 multi-component: every rename applied" begin
        data = load_fixture("sample_v0.0.2_multi.json")
        migrate!(data)
        assert_current(data)

        # top-level keys renamed
        for old in ("Sample", "Buffer", "NMR Tube", "Laboratory Reference", "Notes")
            @test !haskey(data, old)
        end
        @test haskey(data, "sample")
        @test haskey(data, "people")
        @test data["people"]["users"] == Any["Alice", "Bob"]

        # per-component renames
        comps = data["sample"]["components"]
        @test length(comps) == 3
        for c in comps
            @test !haskey(c, "Name")
            @test !haskey(c, "Concentration")
            @test !haskey(c, "Unit")
            @test !haskey(c, "Isotopic labelling")
            @test !haskey(c, "Custom labelling")
            @test haskey(c, "name")
            @test haskey(c, "concentration_or_amount")
            @test haskey(c, "unit")
            @test haskey(c, "isotopic_labelling")
            @test haskey(c, "custom_labelling")
            @test haskey(c, "molecular_weight")
            @test haskey(c, "type")
            @test c["type"] === nothing
        end
        @test !any(c["unit"] == "equiv" for c in comps)

        # sample-level renames
        @test !haskey(data["sample"], "Label")
        @test data["sample"]["label"] == "Test v0.0.2 with multiple components"
        @test data["sample"]["physical_form"] == ""  # added in 0.0.3 → 0.1.0

        # buffer field renames
        @test !haskey(data["buffer"], "pH")
        @test data["buffer"]["ph"] == 7.4
        @test !haskey(data["buffer"], "Solvent")
        @test data["buffer"]["solvent"] == "10% D2O"
        @test !haskey(data["buffer"], "Chemical shift reference")
        @test data["buffer"]["chemical_shift_reference"] == "DSS"
        @test !haskey(data["buffer"], "Reference concentration")
        @test data["buffer"]["reference_concentration"] == 10
        @test !haskey(data["buffer"], "Reference unit")
        @test data["buffer"]["reference_unit"] == "uM"
        @test !haskey(data["buffer"], "Custom solvent")
        @test data["buffer"]["custom_solvent"] == ""

        # buffer component renames
        buf_comps = data["buffer"]["components"]
        @test length(buf_comps) == 2
        for bc in buf_comps
            @test !haskey(bc, "Concentration")
            @test !haskey(bc, "Unit")
            @test haskey(bc, "concentration")
            @test haskey(bc, "unit")
        end

        # NMR tube renames
        @test !haskey(data["nmr_tube"], "Type")
        @test data["nmr_tube"]["type"] == "shigemi"
        @test !haskey(data["nmr_tube"], "Sample Volume (μL)")
        @test data["nmr_tube"]["sample_volume_uL"] == 300
        @test !haskey(data["nmr_tube"], "samplejet_rack_position")  # removed in 0.0.3 → 0.1.0
        @test !haskey(data["nmr_tube"], "samplejet_rack_id")
        @test data["nmr_tube"]["rack_id"] == "rack-001"

        # diameter mapped from "5 mm" to 5.0 and renamed to diameter_mm
        @test !haskey(data["nmr_tube"], "diameter")
        @test data["nmr_tube"]["diameter_mm"] == 5.0

        # reference field renames
        @test !haskey(data["reference"], "Labbook Entry")
        @test data["reference"]["labbook_entry"] == "page 42"
        @test !haskey(data["reference"], "Experiment ID")
        @test data["reference"]["sample_id"] == "EXP-2024-001"
    end

    @testset "v0.0.3 diameter string-to-number mappings" begin
        for (diameter_str, expected) in [("1.7 mm", 1.7), ("3 mm", 3.0), ("5 mm", 5.0), ("", nothing)]
            data = Dict{String,Any}(
                "sample" => Dict{String,Any}(),
                "nmr_tube" => Dict{String,Any}("diameter" => diameter_str),
                "metadata" => Dict{String,Any}("schema_version" => "0.0.3"),
            )
            migrate!(data)
            @test !haskey(data["nmr_tube"], "diameter")
            @test data["nmr_tube"]["diameter_mm"] === expected
        end
    end

    @testset "v0.0.3 samplejet rename, remove, and physical_form addition" begin
        data = Dict{String,Any}(
            "sample" => Dict{String,Any}(),
            "nmr_tube" => Dict{String,Any}(
                "diameter" => "5 mm",
                "samplejet_rack_id" => "rack-001",
                "samplejet_rack_position" => "A3",
            ),
            "metadata" => Dict{String,Any}("schema_version" => "0.0.3"),
        )
        migrate!(data)
        @test data["nmr_tube"]["rack_id"] == "rack-001"
        @test !haskey(data["nmr_tube"], "samplejet_rack_id")
        @test !haskey(data["nmr_tube"], "samplejet_rack_position")
        @test data["sample"]["physical_form"] == ""
    end

    @testset "v0.1.0 chain migrates to current" begin
        data = Dict{String,Any}(
            "sample" => Dict{String,Any}(
                "physical_form" => "solution",
                "components" => Any[
                    Dict{String,Any}(
                    "concentration_or_amount" => 1.0,
                    "unit" => "mM",
                    "isotopic_labelling" => "13C,15N",
                ),
                ],
            ),
            "nmr_tube" => Dict{String,Any}("diameter" => 5.0),
            "metadata" => Dict{String,Any}("schema_version" => "0.1.0"),
        )
        migrate!(data)
        assert_current(data)
        @test !haskey(data["nmr_tube"], "diameter")
        @test data["nmr_tube"]["diameter_mm"] == 5.0
    end

    @testset "v0.2.0 multi-component: wildcard map + set over arrays" begin
        data = load_fixture("sample_v0.2.0_multi.json")
        migrate!(data)
        assert_current(data)

        comps = data["sample"]["components"]
        @test length(comps) == 3
        for c in comps
            @test c["unit"] != "equiv"
            @test haskey(c, "molecular_weight")
            @test haskey(c, "type")
            @test c["type"] === nothing
        end

        # nmr_tube.diameter renamed to diameter_mm in 0.2.0 → 0.3.0
        @test !haskey(data["nmr_tube"], "diameter")
        @test haskey(data["nmr_tube"], "diameter_mm")
        @test data["nmr_tube"]["diameter_mm"] == 5.0
    end

    @testset "v0.3.0 multi-component: type added to every component" begin
        data = load_fixture("sample_v0.3.0_multi.json")
        migrate!(data)
        assert_current(data)

        comps = data["sample"]["components"]
        @test length(comps) == 3
        for c in comps
            @test haskey(c, "type")
            @test c["type"] === nothing
            @test haskey(c, "molecular_weight")
            @test haskey(c, "isotopic_labelling")
        end
    end

    @testset "Edge cases: empty and missing components" begin
        data = load_fixture("sample_v0.2.0_empty_components.json")
        migrate!(data)
        assert_current(data)
        @test data["sample"]["components"] == Any[]

        data = load_fixture("sample_v0.2.0_no_components.json")
        migrate!(data)
        assert_current(data)
        @test !haskey(data["sample"], "components")
    end

    @testset "Already current is no-op and idempotent" begin
        data = load_fixture("sample_v0.4.0_already_current.json")
        before = deepcopy(data)
        migrate!(data)
        @test data == before
        @test data["sample"]["components"][1]["isotopic_labelling"] == "19F"

        # idempotency from an older version
        data = load_fixture("sample_v0.0.2_multi.json")
        migrate!(data)
        first = deepcopy(data)
        migrate!(data)
        @test data == first
    end

    @testset "Stress: 50-element component array" begin
        data = Dict{String,Any}(
            "sample" => Dict{String,Any}(
                "components" => Any[
                    Dict{String,Any}(
                        "Name" => "c$i",
                        "Concentration" => i,
                        "Unit" => iseven(i) ? "mM" : "equiv",
                    ) for i in 1:50
                ],
            ),
            "metadata" => Dict{String,Any}("schema_version" => "0.0.2"),
        )
        migrate!(data)
        comps = data["sample"]["components"]
        @test length(comps) == 50
        for c in comps
            @test haskey(c, "name")
            @test haskey(c, "concentration_or_amount")
            @test c["unit"] != "equiv"
            @test c["type"] === nothing
            @test haskey(c, "molecular_weight")
        end
    end

    @testset "loadsample convenience function" begin
        data = loadsample(joinpath(FIXTURES_DIR, "sample_v0.3.0_multi.json"), PATCH_PATH)
        assert_current(data)
        @test length(data["sample"]["components"]) == 3
    end
end
