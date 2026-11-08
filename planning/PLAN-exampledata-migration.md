# Moving the example data to NMRExampleData.jl

NMRTools 0.3.2 finds its example data with NMRExampleData.jl's `example(name)`, loads them with
`loadnmr` and `loadfid`, and no longer has an `Artifacts.toml`. Its tests and documentation still use the old names (`"2D_HN"` and so on),
which must therefore exist in NMRExampleData as **aliases** of the new datasets. An alias is an
exact match, so these names keep meaning the same data however many datasets are added. Once the
data are in place, the documentation can move to the new names in one sweep.

This branch cannot pass its tests until NMRExampleData has been released and registered with
these datasets. Merge it after that.

## What each dataset must provide

`example(name)` returns the dataset's `path`, so `path` must point at the same thing the old
artifact path did: the Bruker experiment, wherever there is one. The nmrPipe and UCSF files
beside it are declared in the `formats` table of `dataset.toml`, so that the tests and
documentation reach them as `example(name, :nmrpipe)` and `example(name, :ucsf)` rather than
by file name. `loadfid(example(name))` finds the raw data in the experiment folder holding
`path`.

| Old name (alias) | Suggested name | `path` must give | Raw data | `formats` |
|:--|:--|:--|:--|:--|
| `1D_1H` | `1d_1h_ubq` | the experiment folder | yes | |
| `1D_19F` | `1d_19f_<sample>` | the experiment folder | yes | `nmrpipe`: its `test.ft1` |
| `1D_19F_titration` | `titration1d_<sample>` | 11 experiment folders, in order | yes | the folders named `1` to `11`, side by side (the `sumexpts` tests copy `1` to `3` from their parent) |
| `2D_HN` | `2d_hn_ubq` | the experiment folder | yes | `nmrpipe`: its `test.ft2`; `ucsf`: its `hmqc.ucsf` |
| `2D_HN_titration` | `titration2d_<sample>` | `["1/test.ft2", …, "11/test.ft2"]` | **no** | a test checks that `loadfid` fails |
| `pseudo2D_T2` | `relaxation1d_<sample>` | the experiment folder | yes | |
| `pseudo2D_XSTE` | `diffusion1d_<sample>` | the experiment folder | yes | `nmrpipe`: its `test.ft1` |
| `pseudo3D_HN_R2` | `relaxation2d_ubq` | the experiment folder | | |
| `pseudo3D_kinetics` | `kinetics2d_<sample>` | the experiment folder | | |
| `3D_HNCA` | `3d_hnca_ubq` | the experiment folder, with `pdata/131` and `pdata/231` | | `nmrpipe`: the old `3D_HNCA_pipe`'s `smile/test%03d.ft3`; `ucsf`: the old `3D_HNCA_ucsf`'s `hnca.ucsf`, the two the same size |

The BMRB shifts of ubiquitin are no longer bundled: the docs call `loadbmrb(17769)`.

The old artifacts could not be inspected while this was prepared, so the layouts above are
inferred from how the tests and documentation use them. Check each against the old artifact
before adding it.

## Order of adding

Where a shortened name matches several datasets, the first added wins. Add `1d_1h_ubq`,
`1d_19f_…`, `2d_hn_ubq` and `3d_hnca_ubq` before any other dataset starting the same way, so that
`"1d_1h"`, `"1d_19f"`, `"2d_hn"` and `"3d_hnca"` find them.

The old `3D_HNCA_pipe` and `3D_HNCA_ucsf` become formats of `3d_hnca_ubq`, and are no longer
names of their own: the tests and documentation use `example("3D_HNCA", :nmrpipe)`. Copy
their files into the dataset beside the experiment, processed over a narrow region if their
size is a problem, and add them to `formats` by hand.

## Then

1. Add the datasets to NMRExampleData and push to `main`. The Release workflow uploads them and
   registers the new version.
2. Run the NMRTools tests and build the documentation against that version, then merge this
   branch and register NMRTools 0.3.2.
3. Move the documentation and tests to the new names, keeping the aliases for users' code.
