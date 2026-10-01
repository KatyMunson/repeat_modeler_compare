#!/usr/bin/env Rscript
# plot_repeat_compare.R -- three comparison plots from repeat_compare's
# summary tables: (1) class composition (primary "shared" arm) with tissue
# + assembly covariates in the subtitle, (2) divergence landscapes faceted
# by species, (3) shared-vs-own concordance dot plot.
#
# Usage: Rscript plot_repeat_compare.R <class_composition.tsv> \
#          <divergence_landscape.tsv> <arm_concordance.tsv> \
#          <assembly_covariates.tsv> <out_dir>
#
# Pass NONE for <arm_concordance.tsv> to skip plot 3 (shared-arm-only
# report, where there is no own arm to compare against).

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

# Genome-scale bp exceed 32-bit ints; read them as doubles, not integer64
# (whose arithmetic silently breaks when bit64 isn't installed).
class_composition <- fread(class_composition_path, integer64 = "double")
divergence_landscape <- fread(divergence_landscape_path, integer64 = "double")
has_concordance <- arm_concordance_path != "NONE"
assembly_covariates <- fread(assembly_covariates_path, integer64 = "double")

# ---------------------------------------------------------------------------
# One class -> color mapping shared by every plot, so a class keeps its color
# whichever plot it's in and whichever classes a plot happens to contain.
# Listed bottom of the stack first. The eight hues are a validated
# categorical palette, assigned in an order where every pair adjacent in the
# stack stays distinct for common color-vision deficiencies. Unknown, Other
# and Low_complexity are deliberately neutral grays, since they carry no
# biological identity. Classes absent from this list fall back to black, and
# the script warns.
# ---------------------------------------------------------------------------
class_colors <- c(
  DNA            = "#2a78d6",  # blue
  RC             = "#eb6834",  # orange
  LINE           = "#1baf7a",  # aqua
  SINE           = "#eda100",  # yellow
  Retroposon     = "#e87ba4",  # magenta
  LTR            = "#008300",  # green
  Satellite      = "#4a3aa7",  # violet
  Simple_repeat  = "#e34948",  # red
  Low_complexity = "#b0afa9",  # light gray
  Other          = "#5f5e5a",  # dark gray
  Unknown        = "#8f8e88"   # mid gray
)

order_classes <- function(dt) {
  unmapped <- setdiff(unique(dt$class), names(class_colors))
  if (length(unmapped) > 0) {
    warning("classes with no assigned color (drawn black): ", paste(unmapped, collapse = ", "))
    class_colors <<- c(class_colors, setNames(rep("#000000", length(unmapped)), unmapped))
  }
  # ggplot stacks the first factor level on top, so reverse: DNA sits on the
  # baseline and the gray bins on top.
  dt[, class := factor(class, levels = rev(names(class_colors)))]
  dt
}

class_fill <- function() scale_fill_manual(values = class_colors, name = "Class")
class_colour <- function() scale_colour_manual(values = class_colors, name = "Class")

# ---------------------------------------------------------------------------
# 1. Class composition, primary "shared" arm, species side by side.
# ---------------------------------------------------------------------------
shared_composition <- order_classes(class_composition[arm == "shared"])

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
  geom_col(position = "stack", colour = "white", linewidth = 0.4, width = 0.6) +
  class_fill() +
  labs(
    title = "Repeat class composition (shared-library arm)",
    subtitle = covariate_subtitle,
    x = "Species",
    y = "% of non-N assembly length"
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
# 2. Divergence landscapes, primary "shared" arm, faceted by species,
#    stacked by class. Overlap-resolved (each base counted once), so heights
#    are genome fractions and a shared y-axis compares species directly.
# ---------------------------------------------------------------------------
shared_landscape <- divergence_landscape[arm == "shared",
  .(pct_non_n = sum(pct_non_n)), by = .(species, class, kimura_bin)]
shared_landscape <- order_classes(shared_landscape)

p2 <- ggplot(shared_landscape, aes(x = kimura_bin, y = pct_non_n, fill = class)) +
  geom_col(position = "stack", width = 1, colour = "white", linewidth = 0.1) +
  facet_wrap(~species) +
  class_fill() +
  labs(
    title = "Divergence (Kimura) landscape by class (shared-library arm)",
    subtitle = "Overlapping alignments resolved: each base counted once. Simple_repeat and Low_complexity have no divergence.",
    x = "Kimura substitution level (%)",
    y = "% of non-N assembly length"
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
if (has_concordance) {
  arm_concordance <- order_classes(fread(arm_concordance_path, integer64 = "double"))

  p3 <- ggplot(arm_concordance, aes(x = pct_non_n_own, y = pct_non_n_shared, color = class)) +
    geom_point(size = 2) +
    class_colour() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~species) +
    labs(
      title = "Shared-library vs own-library concordance",
      x = "% non-N (own-library arm)",
      y = "% non-N (shared-library arm)"
    ) +
    theme_minimal()

  ggsave(
    file.path(out_dir, "arm_concordance.png"),
    plot = p3,
    width = 8,
    height = 6,
    device = grDevices::png
  )
}

if (interactive()) {
  print(p1)
  print(p2)
  if (has_concordance) print(p3)
}
