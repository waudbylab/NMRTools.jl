# Moving the example data to NMRExampleData.jl

NMRTools 0.3.2 finds its example data with NMRExampleData.jl's `example(name)`, loads them with
`loadnmr` and `loadfid`, and no longer has an `Artifacts.toml`. Its tests and documentation use
the datasets' new names.

This branch cannot pass its tests until NMRExampleData has been released and registered with
these datasets. Merge it after that.

## What each dataset must provide

`example(name)` returns the dataset's `path`, so `path` must point at the same thing the old
artifact path did: the Bruker experiment, wherever there is one. The nmrPipe and UCSF files
beside it are declared in the `formats` table of `dataset.toml`, so that the tests and
documentation reach them as `example(name, :nmrpipe)` and `example(name, :ucsf)` rather than
by file name. `loadfid(example(name))` finds the raw data in the experiment folder holding
`path`.

| Old name | Dataset | `path` | Raw data | `formats` |
|:--|:--|:--|:--|:--|
| `1D_1H` | `1d_1h_ubq` (added) | the experiment folder | yes | |
| `1D_19F` | `1d_19f_fragments` (added; a different spectrum) | the experiment folder | yes | none, so no 1D nmrPipe test |
| `1D_19F_titration` | `1d_19f_titration_sh3` (added) | 11 experiment folders, in order | yes | |
| `2D_HN` | `2d_hn_ubq` (added) | the experiment folder | yes | `nmrpipe`: its `test.ft2`; `ucsf`: its `hmqc.ucsf` |
| `2D_HN_titration` | `titration2d_fir_nbox` (added) | 11 experiment folders, in order | yes | |
| `pseudo2D_T2` | `relaxation1d_aas` (added) | the experiment folder | yes | |
| `pseudo2D_XSTE` | `diffusion1d_asyn` (added) | the experiment folder | yes | `nmrpipe`: its `test.ft1` |
| `pseudo3D_HN_R2` | `relaxation2d_asyn` (added) | the experiment folder | no | |
| `pseudo3D_kinetics` | `kinetics2d_cjun` (added) | the experiment folder | | |
| `3D_HNCA` | `3d_hnca_ubq` (added) | the experiment folder, with `pdata/131` and `pdata/231` | | none, so no 3D nmrPipe or UCSF test |

Every old name is now replaced in the tests and documentation.

The example datasets keep only the real parts of processed spectra, so `allcomponents=true` is
tested on the bundled `test/test-data` (`samples/2`, `15n-trosy-ts3`, `15n-xste-ts3`), and the
docs make a complex spectrum with `autoprocess(...; real=false)`.

The BMRB shifts of ubiquitin are no longer bundled: the docs call `loadbmrb(17769)`.

The old artifacts could not be inspected while this was prepared, so the layouts above are
inferred from how the tests and documentation use them. Check each against the old artifact
before adding it.

## Order of adding

Where a shortened name matches several datasets, the first added wins. Add `1d_1h_ubq`,
`1d_19f_…`, `2d_hn_ubq` and `3d_hnca_ubq` before any other dataset starting the same way, so that
`"1d_1h"`, `"1d_19f"`, `"2d_hn"` and `"3d_hnca"` find them.

The old `3D_HNCA_pipe` and `3D_HNCA_ucsf` become formats of `3d_hnca_ubq`, and are no longer
names of their own: the tests and documentation use `example("3d_hnca_ubq", :nmrpipe)` once it is added. Copy
their files into the dataset beside the experiment, processed over a narrow region if their
size is a problem, and add them to `formats` by hand.

## Then

1. Add the datasets to NMRExampleData and push to `main`. The Release workflow uploads them and
   registers the new version.
2. Move the documentation and tests to the names of the datasets still to be added.
3. Run the NMRTools tests and build the documentation against that version, then merge this
   branch and register NMRTools.
