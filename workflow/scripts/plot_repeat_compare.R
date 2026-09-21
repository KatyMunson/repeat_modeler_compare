#!/usr/bin/env Rscript
# plot_repeat_compare.R -- three comparison plots from repeat_compare's
# summary tables: (1) class composition (primary "shared" arm) with tissue
# + assembly covariates in the subtitle, (2) divergence landscapes faceted
# by species, (3) shared-vs-own concordance dot plot.
#
# Usage: Rscript plot_repeat_compare.R <class_composition.tsv> \
#          <divergence_landscape.tsv> <arm_concordance.tsv> \
#          <assembly_covariates.tsv> <out_dir>

suppressMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
class_composition_path <- args[1]
divergence_landscape_path <- args[2]
arm_concordance_path <- args[3]
assembly_covariates_path <- args[4]
out_dir <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

class_composition <- fread(class_composition_path)
divergence_landscape <- fread(divergence_landscape_path)
arm_concordance <- fread(arm_concordance_path)
assembly_covariates <- fread(assembly_covariates_path)

# ---------------------------------------------------------------------------
# 1. Class composition, primary "shared" arm, species side by side.
# ---------------------------------------------------------------------------
shared_composition <- class_composition[arm == "shared"]

covariate_subtitle <- paste(
  sprintf(
    "%s (%s): N50=%s, %s contigs, %.2f Gb non-N",
    assembly_covariates$species_id,
    assembly_covariates$tissue,
    format(assembly_covariates$contig_n50, big.mark = ","),
    format(assembly_covariates$contig_count, big.mark = ","),
    assembly_covariates$non_n_bp / 1e9
  ),
  collapse = "  |  "
)

p1 <- ggplot(shared_composition, aes(x = species, y = pct_non_n, fill = class)) +
  geom_col(position = "stack") +
  labs(
    title = "Repeat class composition (shared-library arm)",
    subtitle = covariate_subtitle,
    x = "Species",
    y = "% of non-N assembly length",
    fill = "Class"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "class_composition_shared.png"),
  plot = p1,
  width = 8,
  height = 6,
  device = grDevices::png
)

# ---------------------------------------------------------------------------
# 2. Divergence landscapes, faceted by species, stacked by class.
# ---------------------------------------------------------------------------
p2 <- ggplot(divergence_landscape, aes(x = kimura_bin, y = bp, fill = class)) +
  geom_col(position = "stack") +
  facet_wrap(~species, scales = "free_y") +
  labs(
    title = "Divergence (Kimura) landscape by class",
    x = "Kimura substitution level (%)",
    y = "bp",
    fill = "Class"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "divergence_landscape.png"),
  plot = p2,
  width = 10,
  height = 6,
  device = grDevices::png
)

# ---------------------------------------------------------------------------
# 3. Shared vs own concordance dot plot.
# ---------------------------------------------------------------------------
p3 <- ggplot(arm_concordance, aes(x = pct_non_n_own, y = pct_non_n_shared, color = class)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~species) +
  labs(
    title = "Shared-library vs own-library concordance",
    x = "% non-N (own-library arm)",
    y = "% non-N (shared-library arm)",
    color = "Class"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "arm_concordance.png"),
  plot = p3,
  width = 8,
  height = 6,
  device = grDevices::png
)

if (interactive()) {
  print(p1)
  print(p2)
  print(p3)
}
