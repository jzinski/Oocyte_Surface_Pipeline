# Oocyte Surface Pipeline

3D image-analysis pipeline for measuring how two RNA species (*cyclinB1*, *dazl*)
are distributed across the surface of the zebrafish oocyte, comparing wild-type,
Δ6, and Δ6+11 conditions.

## What it does

Given 3-channel confocal z-stacks (C1 = oocyte outline / DAPI, C2 = cyclinB1, C3 = dazl)
the pipeline:

1. Segments the oocyte from C1 in 3D (Gaussian smooth → binarize → keep largest
   connected component → per-Z convex-hull refinement → permuted-axis pass to
   patch slices missed under coarse Z sampling).
2. Sets per-imaging-date C2/C3 intensity thresholds (median + 3·std of
   in-oocyte voxel intensities) and builds C2/C3 masks.
3. Rotates each oocyte so the cyclinB1 centroid defines the animal pole,
   projects into spherical coordinates, and bins mask voxels into 5° angle-from-
   animal-pole elevation shells, producing per-angle surface-coverage curves
   (`elePercC2`, `elePercC3`).
4. Fits a bounded logistic to each coverage curve to extract a half-max
   boundary angle per oocyte, and plots mean ± SD curves and boundary-angle
   violin plots by genotype.

## Folder layout

```
Oocyte_Surface_Pipeline/
├── scripts/                         <- run these, in order
│   ├── step1_generate_meta_files.m  <- C1 segmentation + intensity stats
│   ├── step2_post_meta_thresh.m     <- C2/C3 thresholds + angle binning
│   ├── step3_plot_script.m          <- MATLAB plots + sigmoid fits
│   └── step4_plot_publication.R     <- R/ggplot publication plots
├── functions/                       <- custom MATLAB functions
├── thirdparty/                      <- 3rd-party MATLAB libs (keep licenses)
│   ├── shadedErrorBar/
│   ├── Violinplot-Matlab/
│   └── gifplayer/
└── data/                            <- put image stacks + metadata xlsx here
    ├── Oocyte_sample_list.xlsx      <- metadata template (edit to match your data)
    └── <your>.tif                   <- one 3-channel z-stack per oocyte
```

## Requirements

- MATLAB R2021a or later, with the **Image Processing Toolbox** and
  **Curve Fitting Toolbox**.
- R ≥ 4.0 (for the publication plots) with packages:
  `readxl, dplyr, tidyr, purrr, stringr, ggplot2, patchwork, R.matlab,
  minpack.lm, forcats, scales`.

## Input data format

- Each oocyte is a single `.tif` 3D z-stack with 3 channels interleaved by slice
  (channel order: C1 oocyte outline, C2 cyclinB1, C3 dazl). If starting from
  `.czi`, convert to multi-channel `.tif` in Fiji/ImageJ first.
- Expected pixel size in `xy` = 0.3321 µm, `z` = 2.6 µm. Change at the top of
  each script if your acquisition differs.
- `data/Oocyte_sample_list.xlsx` must contain (at least) these columns:
  `name` (file base, no `.tif`), `date` (integer / yyyymmdd), `genotype`
  (e.g. `wt`, `d6`, `d6+11`).

## Running the pipeline

1. Place all `.tif` stacks and the metadata xlsx in `data/`.
2. In MATLAB, open and run `scripts/step1_generate_meta_files.m`. Output:
   `*_c1mask.mat` per oocyte + `Oocyte_sample_list_thresh.xlsx`.
3. Run `scripts/step2_post_meta_thresh.m`. Output: `*_out.mat`, `*_c2.gif`,
   `*_c3.gif`, `*_masterplot_all.tif` per oocyte, plus
   `Oocyte_sample_list_thresh2.xlsx`.
4. Run `scripts/step3_plot_script.m` for exploratory plots + sigmoid fits.
5. (Optional) In R, open `scripts/step4_plot_publication.R`, `setwd()` to the
   `data/` folder, and source the file for publication-style figures.

Each MATLAB script auto-adds `../functions` and `../thirdparty` to the path and
`cd`s into `../data` on launch, so you can run them from anywhere.

## Output files (per oocyte, in `data/`)

| File | Meaning |
| --- | --- |
| `*_c1mask.mat` | 3D logical mask of the oocyte (channel 1 segmentation) |
| `*_out.mat` | Per-oocyte metrics: `elePercC2`, `elePercC3`, voxel counts, thresholds |
| `*_c2.gif`, `*_c3.gif` | QC overlays of C2 / C3 masks on the raw stack |
| `*_masterplot_all.tif` | Summary figure per oocyte |
| `Oocyte_sample_list_thresh2.xlsx` | Master table with volumes, intensities, thresholds, fitted boundary angles |

## Third-party code

`thirdparty/` contains unmodified copies of:

- `shadedErrorBar` by Rob Campbell — LGPL v3
- `Violinplot-Matlab` by Bastian Bechtold — BSD 3-clause
- `gifplayer` by Vihang Patil — BSD 2-clause

Their original license files are retained in each subfolder.

## License

Custom code in `scripts/` and `functions/` is released under the MIT License
(add a `LICENSE` file if distributing).
