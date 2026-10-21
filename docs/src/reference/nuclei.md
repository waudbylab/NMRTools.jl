# Nuclei, coherences and Xi ratios

## Nuclei

NMRTools defines the nuclei below as instances of [`Nucleus`](@ref). A nucleus can also be looked
up from a string with [`nucleus`](@ref), e.g. `nucleus("15N")`. The columns give the values
returned by each function for nucleus `n`.

```@eval
using NMRTools, Markdown
fmt(x) = isnothing(x) ? "–" : string(x)
rows = ["| `$n` | $(fmt(spin(n))) | $(fmt(gyromagneticratio(n))) | $(fmt(xi_ratio(n; aqueous=true))) | $(fmt(xi_ratio(n; aqueous=false))) |"
        for n in instances(Nucleus)]
Markdown.parse("""
| Nucleus | Spin: `spin(n)` | γ (rad s⁻¹ T⁻¹): `gyromagneticratio(n)` | Xi ratio, DSS (aqueous): `xi_ratio(n)` | Xi ratio, TMS (organic): `xi_ratio(n; aqueous=false)` |
|:--------|:-----|:----------------|:------------------------|:------------------------|
$(join(rows, "\n"))
""")
```

Xi ratios are the IUPAC ratios of each nucleus's reference frequency to that of ¹H in the
reference compound, used for indirect referencing (see
[Reference chemical shifts](../howto/referencing.md)). DSS is the reference for aqueous samples
and TMS for organic solvents.

## Coherences

Coherences are described by [`SQ`](@ref) (single quantum, on one nucleus) and [`MQ`](@ref)
(multiple quantum, a combination of nuclei with coherence orders), both subtypes of
[`Coherence`](@ref). [`coherenceorder`](@ref) gives the total order:

```@setup nuclei
using NMRTools
```

```@repl nuclei
coherenceorder(SQ(H1))
coherenceorder(MQ(((H1, 1), (C13, 1))))     # double quantum
coherenceorder(MQ(((H1, 1), (C13, -1))))    # zero quantum
```

[`gyromagneticratio`](@ref) of a coherence is the sum over its nuclei, weighted by their orders:

```@repl nuclei
gyromagneticratio(MQ(((H1, 1), (C13, 1))))
```
