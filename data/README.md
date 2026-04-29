# data/

Put one 3-channel `.tif` z-stack per oocyte here. The metadata spreadsheet
`Oocyte_sample_list.xlsx` is already in this folder.

For a quick end-to-end test, download the three sample stacks from Zenodo
([doi:10.5281/zenodo.19895264](https://doi.org/10.5281/zenodo.19895264)) and
drop them into this folder — the spreadsheet is pre-filled to match.

Required columns in the spreadsheet:

- `name` - base filename of the stack, without the `.tif` extension
  (e.g. `WT_dazl_cyclinB1_10` for a file `WT_dazl_cyclinB1_10.tif`)
- `date` - imaging date as a plain integer (the pipeline groups per-date for
  intensity thresholding)
- `genotype` - one of `wt`, `d6`, `d6+11` (or whatever labels you use;
  the plot scripts colour-code by these)

All intermediate and output files (`*_c1mask.mat`, `*_out.mat`, `*_c2.gif`,
`*_c3.gif`, `*_masterplot_all.tif`, and updated `Oocyte_sample_list_thresh*.xlsx`
files) will be written here alongside the input data.
