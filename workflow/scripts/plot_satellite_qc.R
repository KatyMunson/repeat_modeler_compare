#!/usr/bin/env Rscript
# plot_satellite_qc.R -- one panel per species: every satellite-library
# motif's genome-wide copy estimate (x, log10) against its monomer length
# (y, log10), from summary/satellite_library_qc.tsv. Colour = copy tier,
# shape = whether the motif passes every QC criterion. Reference lines:
# min_copies and major_min_copies (vertical), the major_min_bp diagonal
# (copies x length = bp), and the monomer length bounds (horizontal).
# Use it to see whether the copy distribution has a natural gap to set
# satellite.qc.min_copies at (same idea as compare_assemblies_satellites
# stage 02's copy_number_diagnostic.png). Motifs with zero hits are drawn
# at x = 0.5. Reads the per-species QC tables so it is available from the
# early satellite_qc target.
#
# Usage: Rscript plot_satellite_qc.R <out.png> <min_copies> <major_min_copies> \
#          <major_min_bp> <min_len> <max_len> <satellite_library_qc.tsv> [...]

suppressMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_png <- args[1]
min_copies <- as.numeric(args[2])
major_min_copies <- as.numeric(args[3])
major_min_bp <- as.numeric(args[4])
min_len <- as.numeric(args[5])
max_len <- as.numeric(args[6])
qc <- rbindlist(lapply(args[7:length(args)], fread), use.names = TRUE)

qc[, x := pmax(est_copies, 0.5)]
qc[, tier := factor(copy_tier, levels = c("major", "minor", "below_floor"),
                    labels = c("major", "minor", "below floor"))]
qc[, qc_result := factor(ifelse(pass_all == "True" | pass_all == TRUE, "passes all QC", "fails a QC check"),
                         levels = c("passes all QC", "fails a QC check"))]
tier_counts <- qc[, .(label = sprintf("%s: %d motifs (%d major, %d pass all)", species, .N,
                                      sum(copy_tier == "major"),
                                      sum(pass_all == "True" | pass_all == TRUE))), by = species]
qc <- merge(qc, tier_counts, by = "species")

# Reference categorical slots 1-3 (validated for all-pairs use in the
# dataviz reference palette): blue, orange, aqua.
tier_colors <- c("major" = "#2a78d6", "minor" = "#eb6834", "below floor" = "#1baf7a")

y_range <- range(c(qc$monomer_len, min_len, max_len), na.rm = TRUE)
diag_line <- data.table(y = 10^seq(log10(y_range[1]), log10(y_range[2]), length.out = 50))
diag_line[, x := major_min_bp / y]

p <- ggplot(qc, aes(x = x, y = monomer_len)) +
  geom_hline(yintercept = c(min_len, max_len), colour = "#8a8984", linewidth = 0.4, linetype = "dotted") +
  geom_vline(xintercept = min_copies, colour = "#52514e", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = major_min_copies, colour = "#52514e", linewidth = 0.5, linetype = "longdash") +
  geom_line(data = diag_line, aes(x = x, y = y), colour = "#52514e", linewidth = 0.4, linetype = "dotdash") +
  geom_point(aes(colour = tier, shape = qc_result), size = 2.6, stroke = 0.9, alpha = 0.9) +
  scale_colour_manual(values = tier_colors, name = "copy tier", drop = FALSE) +
  scale_shape_manual(values = c("passes all QC" = 16, "fails a QC check" = 1), name = NULL, drop = FALSE) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  facet_wrap(~label, ncol = 1) +
  labs(
    x = "estimated genome-wide copies (satellite-screen bp / monomer length, log scale)",
    y = "monomer length (bp, log scale)",
    title = "Satellite library: copy number vs monomer length",
    subtitle = sprintf(paste0(
      "dashed = min_copies (%s); long-dash = major tier (%s copies); dot-dash = major tier by bp (%s bp);\n",
      "dotted = monomer length bounds (%s-%s bp). Zero-hit motifs drawn at 0.5."),
      format(min_copies, big.mark = ","), format(major_min_copies, big.mark = ","),
      format(major_min_bp, big.mark = ","), min_len, max_len)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e8e7e2", linewidth = 0.3),
    plot.background = element_rect(fill = "#fcfcfb", colour = NA),
    strip.text = element_text(hjust = 0, face = "bold"),
    legend.position = "top",
    plot.subtitle = element_text(colour = "#52514e", size = 9)
  )

n_species <- length(unique(qc$species))
ggsave(out_png, p, width = 9, height = 3.2 + 3.3 * n_species, dpi = 150)
