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
    ├── Oocyte_sample_list.xlsx      <- pre-filled for the 3 Zenodo sample stacks
    └── <your>.tif                   <- one 3-channel z-stack per oocyte
```

## Sample data

Three representative `.tif` z-stacks (one per genotype) are deposited on Zenodo
and not included in the repo because of file size:

[**doi:10.5281/zenodo.19895264**](https://doi.org/10.5281/zenodo.19895264)

Files:

| Genotype | Filename | Size |
| --- | --- | --- |
| WT | `WT_dazl_cyclinB1_10.tif` | 189 MB |
| Δ6 | `d6_dazl_cyclinB1_14.tif` | 139 MB |
| Δ6+11 | `d6+11_dazl_cyclinB1_19.tif` | 151 MB |

Download all three and drop them into `data/`. The included
`Oocyte_sample_list.xlsx` is already populated for these three samples, so the
pipeline runs end-to-end with no editing.

## System requirements

### Hardware
- A standard desktop/laptop is sufficient. ~8 GB RAM recommended (each oocyte
  is a ~150–300 MB 3-channel z-stack held in memory while it is processed).
- No GPU and no non-standard hardware required.

### Software
- **MATLAB** R2021a or later, with two toolboxes:
  - [Image Processing Toolbox](https://www.mathworks.com/products/image.html)
  - [Curve Fitting Toolbox](https://www.mathworks.com/products/curvefitting.html)
  - Runs unmodified on [MATLAB Online](https://matlab.mathworks.com/), so no
    local install is required.
- **R** ≥ 4.0 (only for the optional publication plots in step 4), with the
  CRAN packages listed below.

The pipeline was developed on Windows 10/11 desktop MATLAB and is designed to
run unmodified on MATLAB Online.

## Installation guide

### 1. Get the code from GitHub
Either clone the repository:

```bash
git clone https://github.com/jzinski/Oocyte_Surface_Pipeline.git
```

or download it as a ZIP from
<https://github.com/jzinski/Oocyte_Surface_Pipeline> (green **Code** button →
**Download ZIP**) and unzip it.

No compilation or build step is required — the code is plain MATLAB `.m` files
and one R script. The MATLAB scripts add their own `functions/` and
`thirdparty/` folders to the path automatically on launch.

**Typical install time: under 1 minute** (the time to clone/unzip the repo).

### 2. MATLAB and toolboxes
Install MATLAB and the two required toolboxes from MathWorks:
<https://www.mathworks.com/products/matlab.html>. The Image Processing Toolbox
and Curve Fitting Toolbox are added from the MATLAB installer or **Add-On
Explorer**. Alternatively, use [MATLAB Online](https://matlab.mathworks.com/)
with no install. Typical toolbox install time is a few minutes on a normal
broadband connection (dominated by download, not by this repo).

### 3. R packages (optional — step 4 only)
The publication plots use these CRAN packages:

```r
install.packages(c("readxl", "dplyr", "tidyr", "purrr", "stringr",
                   "ggplot2", "patchwork", "R.matlab", "minpack.lm",
                   "forcats", "scales"))
```

CRAN: <https://cran.r-project.org/>. Typical install time is 1–3 minutes
depending on which packages are already present.

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

## Demo

A small test dataset of three oocytes (one per genotype) is provided on Zenodo
so the pipeline can be run end to end without the full image set.

### Demo data
Download the three `.tif` stacks from
[doi:10.5281/zenodo.19895264](https://doi.org/10.5281/zenodo.19895264) into the
`data/` folder. The bundled `data/Oocyte_sample_list.xlsx` is already filled in
for exactly these three files, so no editing is needed.

| Genotype | File |
| --- | --- |
| WT | `WT_dazl_cyclinB1_10.tif` |
| Δ6 | `d6_dazl_cyclinB1_14.tif` |
| Δ6+11 | `d6+11_dazl_cyclinB1_19.tif` |

### Run the demo
In MATLAB, run the scripts in order:

```matlab
run scripts/step1_generate_meta_files.m   % C1 segmentation + intensity stats
run scripts/step2_post_meta_thresh.m      % C2/C3 thresholds + angle binning
run scripts/step3_plot_script.m           % plots + sigmoid fits
```

To capture the run time, wrap a step with `tic`/`toc`, e.g.:

```matlab
tic; run scripts/step1_generate_meta_files.m; toc
```

### Expected output
For each of the three oocytes the demo writes, into `data/`:

- `*_c1mask.mat` — 3D oocyte mask (from step 1)
- `*_out.mat` — per-oocyte metrics (`elePercC2`, `elePercC3`, voxel counts, thresholds)
- `*_c2.gif`, `*_c3.gif`, `*_c1.gif` — QC overlays of the masks on the raw stack
- `*_masterplot_all.tif` — per-oocyte summary figure

plus the updated tables `Oocyte_sample_list_thresh.xlsx` and
`Oocyte_sample_list_thresh2.xlsx`, and the genotype coverage / half-max plots
from step 3.

The key per-oocyte readout is the half-max angle from the animal pole
(0° = animal pole, 180° = vegetal), reported in the output tables. The expected
biological pattern is that cyclinB1 caps the animal pole (small half-max) while
dazl fills the vegetal (large half-max), with the cyclinB1 cap broadest in
Δ6+11. The `*_masterplot_all.tif` figure for each oocyte is the easiest way to
confirm a run at a glance (raw image, coverage heatmaps, fraction-coverage
curves).

> The exact expected values and run time for this demo are being captured from
> a clean run on MATLAB Online and will be filled in here. *(to complete)*

### Expected run time
*To be filled in from the demo run on MATLAB Online* (wrap each step in
`tic`/`toc` as shown above and record the total for the three oocytes).

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
