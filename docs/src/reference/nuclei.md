# Nuclei, coherences and Xi ratios

## Nuclei

NMRTools defines the nuclei below as instances of [`Nucleus`](@ref). Their spin quantum numbers,
gyromagnetic ratios and Xi ratios are returned by [`spin`](@ref), [`gyromagneticratio`](@ref) and
[`xi_ratio`](@ref). A nucleus can also be looked up from a string with [`nucleus`](@ref), e.g.
`nucleus("15N")`.

```@example nuclei
using NMRTools, Markdown
fmt(x) = isnothing(x) ? "–" : string(x)
rows = ["| `$n` | $(fmt(spin(n))) | $(fmt(gyromagneticratio(n))) | $(fmt(xi_ratio(n; aqueous=true))) | $(fmt(xi_ratio(n; aqueous=false))) |"
        for n in instances(Nucleus)]
Markdown.parse("""
| Nucleus | Spin | γ (rad s⁻¹ T⁻¹) | Xi ratio (DSS, aqueous) | Xi ratio (TMS, organic) |
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

```@example nuclei
coherenceorder(SQ(H1)), coherenceorder(MQ(((H1, 1), (C13, 1)))), coherenceorder(MQ(((H1, 1), (C13, -1))))
```
