# ─────────────────────────────────────────────────────────────────────────────
# Oocyte plots in R (no saving)
# Replicates the MATLAB plotting pipeline with tidyverse + ggplot2
# Requires: readxl, dplyr, tidyr, purrr, stringr, ggplot2, patchwork, R.matlab, minpack.lm
#
# Before running: set the working directory to ../data (the folder containing the
# _out.mat files and the Oocyte_sample_list_thresh2.xlsx produced by step 2).
#   setwd("../data")        # from this scripts/ folder, or
#   setwd("<absolute path to data>")
# ─────────────────────────────────────────────────────────────────────────────

# 0) Libraries -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(R.matlab)     # readMat()
  library(minpack.lm)   # nlsLM() with bounds
})

# 1) Config (edit as needed) ---------------------------------------------------
tableName   <- "Oocyte_sample_list_thresh2"

# pixel sizes (unused in plots but kept here for parity with MATLAB)
xsize <- 0.3321
ysize <- 0.3321
zsize <- 2.6

qualThresh <- 0

# colors from MATLAB (convert to hex)
genoColorSet <- rbind(c(.75,0,0), c(0,.75,0), c(0,0,.75))
dateColorSet <- rbind(c(.75,.75,.75), c(.5,.5,.5), c(.25,.25,.25), c(.75,0,.75))

to_hex <- function(mat) apply(mat, 1, function(x) rgb(x[1], x[2], x[3], 1))
geno_cols <- to_hex(genoColorSet)
date_cols <- to_hex(dateColorSet)

angles <- seq(0, 175, by = 5)
angle_mid <- angles + 2.5

# 2) Read metadata and prune ---------------------------------------------------
metaAll <- read_xlsx(paste0(tableName, ".xlsx")) %>%
  # keep column names as-is to match MATLAB (e.g., good_C1_ID, C3_Thresh, name, genotype, date)
  mutate(
    genotype = as.character(genotype),
    name     = as.character(name)
  )

metaAllPrune <- metaAll %>%
  filter(.data$good_C1_ID > qualThresh, .data$C3_Thresh != 0)

stopifnot(nrow(metaAllPrune) > 0)

# 3) Ingest *_out.mat for each sample -----------------------------------------
# Each .mat must contain: elePercC2, elePercC3, H1counts, H1countsC2, Cencounts
read_one <- function(base) {
  mat_path <- paste0(base, "_out.mat")
  if (!file.exists(mat_path)) stop("Missing file: ", mat_path)
  m <- readMat(mat_path)
  # Coerce to numeric row vectors; pad/trim to match 'angles'
  as_num <- function(x) {
    if (is.null(x)) return(rep(NA_real_, length(angles)))
    v <- as.numeric(x)
    if (length(v) == length(angles)) return(v)
    if (length(v) < length(angles)) return(c(v, rep(NA_real_, length(angles) - length(v))))
    v[seq_along(angles)]
  }
  list(
    elePercC2   = as_num(m$elePercC2),
    elePercC3   = as_num(m$elePercC3),
    H1counts    = as_num(m$H1counts),
    H1countsC2  = as_num(m$H1countsC2),
    Cencounts   = as_num(m$Cencounts)
  )
}

mat_list <- map(metaAllPrune$name, read_one)

# 4) Build long data frame of per-angle values --------------------------------
df_long <- map2_dfr(mat_list, seq_len(nrow(metaAllPrune)), function(m, i) {
  tibble(
    name     = metaAllPrune$name[i],
    genotype = metaAllPrune$genotype[i],
    date     = metaAllPrune$date[i],
    angle    = angles,
    elePercC2 = m$elePercC2,
    elePercC3 = m$elePercC3,
    H1counts  = m$H1counts,
    H1countsC2= m$H1countsC2,
    Cencounts = m$Cencounts
  )
}) %>%
  mutate(
    H1volume_all  = H1counts  / Cencounts,
    H1volumeC2_all= H1countsC2 / Cencounts
  )

# Helpers for mean + sd ribbons
summarise_curve <- function(df, y) {
  df %>%
    group_by(genotype, angle, date) %>%
    summarise(
      mean = mean(.data[[y]], na.rm = TRUE),
      sd   = sd(.data[[y]], na.rm = TRUE),
      .groups = "drop"
    )
}

# 5) Per-date 2×2 figure (C2, C3, H1volumeC2_all, H1volume_all) ---------------
uniqueDates <- sort(unique(df_long$date))
uniqueGeno  <- unique(df_long$genotype)
if (length(geno_cols) < length(uniqueGeno)) {
  # extend gently if >3 genotypes
  geno_cols <- scales::hue_pal()(length(uniqueGeno))
}

for (d in uniqueDates) {
  d_sub <- df_long %>% filter(date == d)
  
  # C2
  c2_sum <- summarise_curve(d_sub, "elePercC2")
  p1 <- ggplot() +
    geom_line(data = d_sub, aes(angle, elePercC2, group = interaction(name, genotype), color = genotype),
              linewidth = 0.6, alpha = 0.35) +
    geom_ribbon(data = c2_sum, aes(angle, ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = genotype),
                alpha = 0.15) +
    geom_line(data = c2_sum, aes(angle, mean, color = genotype), linewidth = 1.5) +
    scale_color_manual(values = geno_cols) +
    scale_fill_manual(values = geno_cols) +
    coord_cartesian(xlim = c(0,180), ylim = c(0,1)) +
    labs(x = "Angle From Animal Pole (degrees)", y = "Fraction Surface Covered by CyclinB1",
         title = "CyclinB1 Fraction Coverage") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  
  # C3
  c3_sum <- summarise_curve(d_sub, "elePercC3")
  p2 <- ggplot() +
    geom_line(data = d_sub, aes(angle, elePercC3, group = interaction(name, genotype), color = genotype),
              linewidth = 0.6, alpha = 0.35) +
    geom_ribbon(data = c3_sum, aes(angle, ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = genotype),
                alpha = 0.15) +
    geom_line(data = c3_sum, aes(angle, mean, color = genotype), linewidth = 1.5) +
    scale_color_manual(values = geno_cols) +
    scale_fill_manual(values = geno_cols) +
    coord_cartesian(xlim = c(0,180), ylim = c(0,1)) +
    labs(x = "Angle From Animal Pole (degrees)", y = "Fraction Surface Covered by Dazl",
         title = "Dazl Fraction Coverage") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  
  # H1volumeC2_all
  vc2_sum <- summarise_curve(d_sub, "H1volumeC2_all")
  p3 <- ggplot() +
    geom_line(data = d_sub, aes(angle, H1volumeC2_all, group = interaction(name, genotype), color = genotype),
              linewidth = 0.6, alpha = 0.35) +
    geom_line(data = vc2_sum, aes(angle, mean, color = genotype), linewidth = 1.5) +
    scale_color_manual(values = geno_cols) +
    labs(x = "Angle From Animal Pole (degrees)", y = "Percent Volume C2") +
    theme_minimal(base_size = 12) + theme(legend.position = "none") +
    coord_cartesian(xlim = c(0,180))
  
  # H1volume_all
  v_sum <- summarise_curve(d_sub, "H1volume_all")
  p4 <- ggplot() +
    geom_line(data = d_sub, aes(angle, H1volume_all, group = interaction(name, genotype), color = genotype),
              linewidth = 0.6, alpha = 0.35) +
    geom_line(data = v_sum, aes(angle, mean, color = genotype), linewidth = 1.5) +
    scale_color_manual(values = geno_cols) +
    labs(x = "Angle From Animal Pole (degrees)", y = "Percent Volume C3") +
    theme_minimal(base_size = 12) + theme(legend.position = "none") +
    coord_cartesian(xlim = c(0,180))
  
  ((p1 + p2) / (p3 + p4)) + plot_annotation(title = paste("Date =", d))
}

# Colors mapped to new labels: WT (green), Δ6 (blue), Δ6+11 (yellow)
geno_cols_named <- c("WT" = "#2E8B57", "Δ6" = "#1F78B4", "Δ6+11" = "#FFD23F")

LINE_W <- 2.2; RIB_A <- 0.25; BASE <- 18

# Recode + order genotypes to WT, Δ6, Δ6+11
keep_dates <- uniqueDates[1:2]
d_sub <- df_long %>%
  dplyr::filter(date %in% keep_dates) %>%
  dplyr::mutate(
    angle_mid = angle + 2.5,
    genotype = forcats::fct_recode(
      factor(genotype),
      "WT"   = "wt",
      "Δ6"   = "d6",
      "Δ6+11"= "d6+11"
    ) %>% forcats::fct_relevel("WT","Δ6","Δ6+11")
  )

theme_pub_base <- theme_minimal(base_size = BASE) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = BASE),
    axis.text  = element_text(size = BASE - 2),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), color = "grey70"),
    legend.title = element_text(size = BASE - 3),
    legend.text  = element_text(size = BASE - 3),
    plot.title = element_blank()
  )

# summaries
c2_sum <- d_sub %>% dplyr::group_by(genotype, angle_mid) %>%
  dplyr::summarise(mean = mean(elePercC2, na.rm = TRUE),
                   sd   = sd(elePercC2,   na.rm = TRUE), .groups = "drop")
c3_sum <- d_sub %>% dplyr::group_by(genotype, angle_mid) %>%
  dplyr::summarise(mean = mean(elePercC3, na.rm = TRUE),
                   sd   = sd(elePercC3,   na.rm = TRUE), .groups = "drop")

# --- CyclinB1 (legend top-right) ---
p_pub_c2 <- ggplot(c2_sum, aes(angle_mid, mean, color = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
              alpha = RIB_A, color = NA) +
  geom_line(linewidth = LINE_W) +
  scale_color_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
  scale_fill_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
  scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  labs(
    x = "Angle From Animal Pole (degrees)",
    y = expression(paste("Fraction Surface Covered by ", italic(cyclinB1)))
  ) +
  guides(color = guide_legend(override.aes = list(fill = NA, linewidth = LINE_W))) +
  theme_pub_base +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

# --- Dazl (legend top-left) ---
p_pub_c3 <- ggplot(c3_sum, aes(angle_mid, mean, color = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
              alpha = RIB_A, color = NA) +
  geom_line(linewidth = LINE_W) +
  scale_color_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
  scale_fill_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
  scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  labs(
    x = "Angle From Animal Pole (degrees)",
    y = expression(paste("Fraction Surface Covered by ", italic(dazl)))
  ) +
  guides(color = guide_legend(override.aes = list(fill = NA, linewidth = LINE_W))) +
  theme_pub_base +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

p_pub_c2
p_pub_c3
# (p_pub_c2 | p_pub_c3)  # side-by-side if you want



# 6) “First two dates for pub”: means + ribbons, genotypes overlaid ------------
if (length(uniqueDates) >= 2) {
  keep_dates <- uniqueDates[1:2]
  d_sub <- df_long %>%
    dplyr::filter(date %in% keep_dates) %>%
    dplyr::mutate(
      angle_mid = angle + 2.5,
      # ensure desired legend order + mapping
      genotype = forcats::fct_recode(
        factor(genotype),
        "wt"    = "wt",
        "d6"    = "d6",
        "d6+11" = "d6+11"
      ) %>% forcats::fct_relevel("wt","d6","d6+11")
    )
  
  # color map: wt = green, d6 = blue, d6+11 = yellow
  geno_cols_named <- c("wt" = "#2E8B57",    # sea green
                       "d6" = "#1F78B4",    # blue
                       "d6+11" = "#FFD23F") # yellow
  
  theme_pub <- theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85"),
      legend.position = c(0.18, 0.80),
      legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                       color = "grey80"),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 11),
      axis.title   = element_text(size = 13),
      axis.text    = element_text(size = 12),
      plot.tag     = element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 0.98)
    )
  
  # summarise mean ± SD
  c2_sum <- d_sub %>%
    dplyr::group_by(genotype, angle_mid) %>%
    dplyr::summarise(mean = mean(elePercC2, na.rm = TRUE),
                     sd   = sd(elePercC2,   na.rm = TRUE), .groups = "drop")
  
  c3_sum <- d_sub %>%
    dplyr::group_by(genotype, angle_mid) %>%
    dplyr::summarise(mean = mean(elePercC3, na.rm = TRUE),
                     sd   = sd(elePercC3,   na.rm = TRUE), .groups = "drop")
  
  p_pub_c2 <- ggplot(c2_sum, aes(angle_mid, mean,
                                 color = genotype, fill = genotype)) +
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                alpha = 0.20, color = NA) +
    geom_line(linewidth = 1.4) +
    scale_color_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
    scale_fill_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
    scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    labs(x = "Angle From Animal Pole (degrees)",
         y = "Fraction Surface Covered by CyclinB1",
         title = "Surface Covered by CyclinB1", tag = "B") +
    guides(color = guide_legend(override.aes = list(fill = NA, alpha = 1, linewidth = 1.4))) +
    theme_pub + theme(plot.title = element_text(color = "#2E8B57", face = "bold", hjust = 0.5))
  
  p_pub_c3 <- ggplot(c3_sum, aes(angle_mid, mean,
                                 color = genotype, fill = genotype)) +
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                alpha = 0.20, color = NA) +
    geom_line(linewidth = 1.4) +
    scale_color_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
    scale_fill_manual(values = geno_cols_named, breaks = names(geno_cols_named)) +
    scale_x_continuous(breaks = seq(0, 180, 20), limits = c(0, 180), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    labs(x = "Angle From Animal Pole (degrees)",
         y = "Fraction Surface Covered by Dazl",
         title = "Surface Covered by Dazl", tag = "B′") +
    guides(color = guide_legend(override.aes = list(fill = NA, alpha = 1, linewidth = 1.4))) +
    theme_pub + theme(plot.title = element_text(color = "#D62728", face = "bold", hjust = 0.5))
  
  # show side-by-side (no saving)
  p_pub_c2
  p_pub_c3
  # If you want them together:
  # (p_pub_c2 | p_pub_c3)
}


# 7) “All dates on same plots”: per-genotype means per date --------------------
# Reuse dateColorSet (will recycle if >4 dates)
p_list <- list()
for (g in uniqueGeno) {
  g_sub <- df_long %>% filter(genotype == g)
  c2m <- g_sub %>%
    group_by(date, angle) %>%
    summarise(mean = mean(elePercC2, na.rm = TRUE), .groups = "drop")
  
  p_c2 <- ggplot(c2m, aes(angle, mean, color = factor(date))) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = rep(date_cols, length.out = length(uniqueDates))) +
    coord_cartesian(xlim = c(0,180)) +
    labs(x = "Angle From Animal Pole (degrees)", y = "Fraction Surface Covered by CyclinB1",
         title = paste("CyclinB1 Fraction Coverage —", g), color = "Date") +
    theme_minimal(base_size = 12)
  
  c3m <- g_sub %>%
    group_by(date, angle) %>%
    summarise(mean = mean(elePercC3, na.rm = TRUE), .groups = "drop")
  
  p_c3 <- ggplot(c3m, aes(angle, mean, color = factor(date))) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = rep(date_cols, length.out = length(uniqueDates))) +
    coord_cartesian(xlim = c(0,180)) +
    labs(x = "Angle From Animal Pole (degrees)", y = "Fraction Surface Covered by Dazl",
         title = paste("Dazl Fraction Coverage —", g), color = "Date") +
    theme_minimal(base_size = 12)
  
  p_list <- append(p_list, list(p_c2, p_c3))
}
wrap_plots(p_list, ncol = 2) + plot_annotation(title = "Date Compare by Genotype")

# 8) Sigmoid fits per sample (C2 increasing, C3 mirrored) ----------------------
# Logistic forms (match MATLAB):
# C2: c3 / (1 + exp(-c1 * (x - c2)))   with c1 in [-1, -0.09], c2 in [5,70], c3 in [0,1]
# C3: c3 / (1 + exp(-c1 * (-x + c2)))  (mirrored in x)
logit_inc <- function(x, c1, c2, c3) c3 / (1 + exp(-c1 * (x - c2)))
logit_mir <- function(x, c1, c2, c3) c3 / (1 + exp(-c1 * (-x + c2)))

fit_one <- function(x, y, mirrored = FALSE) {
  df <- tibble(x = x, y = y) %>% filter(is.finite(x), is.finite(y))
  if (nrow(df) < 8) return(NA)  # too few points
  start <- list(c1 = -0.5, c2 = if (mirrored) 120 else 30, c3 = 1)
  lower <- c(c1 = -1,   c2 = 5,   c3 = 0)
  upper <- c(c1 = -0.09, c2 = 170, c3 = 1)
  f <- tryCatch(
    if (!mirrored) {
      nlsLM(y ~ logit_inc(x, c1, c2, c3), data = df,
            start = start, lower = lower, upper = upper, control = nls.lm.control(maxiter = 200))
    } else {
      nlsLM(y ~ logit_mir(x, c1, c2, c3), data = df,
            start = start, lower = lower, upper = upper, control = nls.lm.control(maxiter = 200))
    },
    error = function(e) NULL
  )
  if (is.null(f)) return(NA)
  coef(f)[["c2"]]  # boundary angle
}

# Compute boundaries per sample
bounds_df <- df_long %>%
  group_by(name, genotype, date) %>%
  summarise(
    C2_bound = fit_one(angles, elePercC2, mirrored = FALSE),
    C3_bound = fit_one(angles, elePercC3, mirrored = TRUE),
    .groups = "drop"
  )

# Quick visual: per (genotype × date) facet with all sample curves + means -----
plot_sig_panel <- function(var, title_txt, ylab_txt) {
  base <- df_long %>%
    left_join(bounds_df, by = c("name","genotype","date"))
  
  means <- base %>%
    group_by(genotype, date, angle) %>%
    summarise(mean = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  ggplot() +
    geom_line(data = base, aes(angle, .data[[var]], group = name, color = genotype), alpha = 0.35, linewidth = 0.6) +
    geom_line(data = means, aes(angle, mean, color = genotype), linewidth = 1.2) +
    scale_color_manual(values = geno_cols) +
    coord_cartesian(xlim = c(0,180), ylim = c(0,1)) +
    labs(x = "Angle From Animal Pole (degrees)", y = ylab_txt, title = title_txt) +
    theme_minimal(base_size = 12) +
    facet_grid(genotype ~ date, scales = "fixed")
}

plot_sig_panel("elePercC2", "Fraction Coverage (C2)", "Fraction Surface Covered by CyclinB1")
plot_sig_panel("elePercC3", "Fraction Coverage (C3)", "Fraction Surface Covered by Dazl")
