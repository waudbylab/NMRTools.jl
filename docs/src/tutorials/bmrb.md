# Assigning ubiquitin from the BMRB

The [BMRB](https://bmrb.io) holds assigned chemical shifts for thousands of proteins. When your
protein has an entry, the shifts predict where each residue's peak should appear in any
correlation experiment. This tutorial uses the shifts of ubiquitin (BMRB entry 17769) to assign an
HNCA spectrum of ubiquitin: predict the ¹H,¹⁵N peaks, register them onto the spectrum, and check
the assignments in the HNCA's ¹³C dimension.

```@example bmrb
using NMRTools, CairoMakie
```

## The shifts

[`loadbmrb`](@ref NMRTools.NMRIO.loadbmrb) downloads the shifts of a BMRB entry, given its number:

```julia
shifts = loadbmrb(17769)
```

The shifts of entry 17769 come with NMRTools as example data, so this page runs without network
access:

```@example bmrb
shifts = exampledata("BMRB_ubiquitin")
shifts[1:5]
```

Each [`ChemicalShift`](@ref) records the residue number and name, the atom, its nucleus and the
shift in ppm. The entry has 1001 shifts:

```@example bmrb
length(shifts), shifts[1].residuename, shifts[1].atom, shifts[1].shift
```

## Predicting the peaks

[`peaklist`](@ref) assembles shifts into the peaks of an experiment. The arguments name the atoms
correlated in each dimension, in spectrum order with the direct dimension first. For the backbone
amide peaks of an ¹H,¹⁵N HSQC:

```@example bmrb
predicted = peaklist(shifts, "H", "N")
predicted[1:3]
```

The labels give the residue and atoms. Plotted on their own, the predicted peaks are a simulated
HSQC:

```@example bmrb
fig, ax, _ = nmrplot(predicted)
save("bmrb-predicted.svg", fig); nothing # hide
```

![](bmrb-predicted.svg)

## Comparing with a spectrum

The HNCA was recorded at 298 K, the temperature of the BMRB entry. Its ¹H,¹⁵N projection, the
maximum over the ¹³C dimension, looks like an HSQC:

```@example bmrb
hnca = exampledata("3D_HNCA")
hn_proj = dropdims(maximum(hnca; dims=F3Dim); dims=F3Dim) / 3
fig, ax, _ = nmrplot(hn_proj; title="HNCA projection with BMRB peaks")
nmrplot!(ax, predicted)
save("bmrb-before.svg", fig); nothing # hide
```

![](bmrb-before.svg)

Contour levels are set from the noise level of the 3D spectrum, which the projection keeps in its
metadata. Taking the maximum over the ¹³C dimension raises the noise floor above that level, so the
projection is divided by 3 to keep the lowest contours above the noise.

The predicted peaks sit close to the observed ones, but are systematically offset: the BMRB shifts
were referenced under the deposited sample conditions, not ours. A few residues, such as G47 and
A46, lie outside the spectral window of the ¹⁵N dimension.

## Registering the prediction

[`register`](@ref) finds the single offset that best matches a peak list
to a spectrum, and returns the list moved by it ([`register!`](@ref) moves the list in place):

```@example bmrb
registered = register(predicted, hn_proj)
fig, ax, _ = nmrplot(hn_proj; title="after registration")
nmrplot!(ax, registered)
save("bmrb-after.svg", fig); nothing # hide
```

![](bmrb-after.svg)

The offset applied is logged. Only a uniform shift is applied: individual peaks are not moved, so
the predicted positions keep their relative pattern.

## Reading off assignments

The registered list is now a set of assignments for the spectrum. Zooming in on a crowded region
shows which residue each peak belongs to (left). [`findpeak`](@ref) snaps each predicted peak onto
the nearest maximum in the spectrum, keeping its label, to give the observed positions (right):

```@example bmrb
region(p) = 7.6 < ppm(p, 1) < 8.4 && 119 < ppm(p, 2) < 124
observed = [findpeak(hn_proj, ppm(p); label=label(p)) for p in filter(region, registered)]
fig = Figure(size=(1000, 450))
ax1 = nmrplot(fig[1, 1], hn_proj[7.6 .. 8.4, 119 .. 124]; title="registered prediction").axis
nmrplot!(ax1, filter(region, registered))
ax2 = nmrplot(fig[1, 2], hn_proj[7.6 .. 8.4, 119 .. 124]; title="observed positions").axis
nmrplot!(ax2, observed)
save("bmrb-zoom.svg", fig); nothing # hide
```

![](bmrb-zoom.svg)

```@example bmrb
observed[1:3]
```

## Side chains

Atom names select specific atoms: `"H", "N"` gives only backbone amides. Passing nuclei instead
gives every directly bonded N–H pair, which is what an HSQC shows, including the side-chain NH₂
groups of Asn and Gln and the side-chain NH of Arg and Trp. Side chains do not appear in the HNCA,
so we compare with the ¹H,¹⁵N spectrum `2D_HN` instead:

```@example bmrb
allnh = peaklist(shifts, H1, N15)
sidechains = filter(p -> !(p in predicted), allnh)
length(predicted), length(allnh)
```

```@example bmrb
nh2(p) = 6.4 < ppm(p, 1) < 7.9 && 108.5 < ppm(p, 2) < 114.5
hsqc = exampledata("2D_HN")
fig, ax, _ = nmrplot(hsqc[6.4 .. 7.9, 108.5 .. 114.5]; title="side-chain NH₂ region")
nmrplot!(ax, filter(nh2, register(sidechains, hsqc)); color=:red)
save("bmrb-sidechains.svg", fig); nothing # hide
```

![](bmrb-sidechains.svg)

## Checking assignments in the HNCA

An HNCA correlates each amide with the Cα of its own residue (`"CA"`) and of the preceding residue
(`"CA-1"`). Predicting both and registering them onto the 3D spectrum gives the expected Cα shifts
for each amide:

```@example bmrb
ca = register(peaklist(shifts, "H", "N", "CA"), hnca)
caprev = register(peaklist(shifts, "H", "N", "CA-1"), hnca)
ca[1:2]
```

A strip is a narrow ¹H range of the HNCA at one amide's ¹⁵N shift. Plotting strips for consecutive
residues, with the predicted Cα shifts marked, shows the sequential walk: each residue's Cα(i−1)
peak lines up with the Cα(i) peak of the strip before it.

```@example bmrb
pick(list, res) = only(filter(p -> startswith(label(p), res * "H-"), list))
residues = ["I13", "T14", "L15", "E16"]
fig = Figure(size=(900, 450))
for (i, res) in enumerate(residues)
    p, q = pick(ca, res), pick(caprev, res)
    h, n = ppm(p, 1), ppm(p, 2)
    strip = hnca[(h - 0.15) .. (h + 0.15), Near(n), :]
    ax = nmrplot(fig[1, i], strip; title=res, xlabel="¹H (ppm)",
                 ylabel=i == 1 ? "¹³C (ppm)" : "").axis
    hlines!(ax, [ppm(p, 3)]; color=:black, linestyle=:dash)
    hlines!(ax, [ppm(q, 3)]; color=:red, linestyle=:dot)
end
save("bmrb-strips.svg", fig); nothing # hide
```

![](bmrb-strips.svg)

Dashed black lines mark the predicted Cα(i), dotted red lines Cα(i−1).

## A spectrum at a different temperature

The ¹H,¹⁵N spectrum `2D_HN` was recorded at 277 K, 21 K below the BMRB entry. Amide shifts change
with temperature by different amounts for each residue, so a uniform offset cannot match every
peak:

```@example bmrb
cold = exampledata("2D_HN")
cold[:temperature]
```

```@example bmrb
moved = register(predicted, cold)
fig, ax, _ = nmrplot(cold; title="2D_HN at 277 K with registered BMRB peaks")
nmrplot!(ax, moved)
save("bmrb-cold.svg", fig); nothing # hide
```

![](bmrb-cold.svg)

Most peaks are close enough to identify, but some need checking by hand or with a 3D experiment.

## Saving the assignments

[`savepeaks`](@ref) writes the registered list as a `peaklist.csv`, the starting point for an
analysis in [NMRAnalysis.jl](https://waudbylab.org/NMRAnalysis.jl):

```julia
savepeaks("ubiquitin-assigned.csv", registered)
```

## Next steps

- [Peak lists from the BMRB](../tools/bmrb.md): every way of writing a correlation, including
  CON, HNCO and methyl spectra
- [Register a peak list onto a spectrum](../howto/register.md)
