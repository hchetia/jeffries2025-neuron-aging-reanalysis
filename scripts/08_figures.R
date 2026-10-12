# Figures for the README. Run from project root:  Rscript scripts/08_figures.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

pal  <- c(housekeeping = "#2a78d6", neuron_specific = "#eb6834")
labs_set <- c(housekeeping = "Housekeeping", neuron_specific = "Neuron-specific")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e0"; surf <- "#fcfcfb"
theme_set(theme_minimal(base_size = 11) + theme(
  plot.background = element_rect(fill = surf, colour = NA), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = grid, linewidth = 0.3), text = element_text(colour = ink),
  axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = 12),
  plot.subtitle = element_text(colour = ink2), legend.position = "top", legend.title = element_blank()))
save <- function(p, f, w, h) ggsave(file.path("results/figures", f), p, width = w, height = h, dpi = 200, bg = surf)

## Fig 1. Expression change vs gene length, excitatory neurons ---------------------
fc <- fread("data/processed/geneset_log2fc_with_length.tsv")[cell_type == "ext"]
win <- fread("results/tables/expression_vs_length_models.tsv")[analysis == "B_length_matched_ext", level][1]
rng <- as.numeric(regmatches(win, gregexpr("[0-9]+", win))[[1]]) * 1e3
p1 <- ggplot(fc, aes(gene_length, log2fc, colour = set)) +
  annotate("rect", xmin = rng[1], xmax = rng[2], ymin = -Inf, ymax = Inf, fill = "#f0efea", colour = NA) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.15, aes(fill = set)) +
  scale_x_log10(labels = function(x) paste0(x / 1e3, " kb")) +
  scale_colour_manual(values = pal, labels = labs_set) + scale_fill_manual(values = pal, guide = "none") +
  labs(x = "Gene length (log scale)", y = "log2(elderly / adult), excitatory neurons",
       title = "Downregulation tracks gene length in both gene classes",
       subtitle = sprintf("Shaded: shared length window (%s). Data: Jeffries et al. 2025, Table S11", sub("window ", "", win)))
save(p1, "fig1_expression_vs_length.png", 7, 5)

## Fig 2. Length slope by cell type -----------------------------------------------
cs <- fread("results/tables/expression_vs_length_models.tsv")[analysis == "C_by_cell_type"]
ord <- c("L2/3", "L4", "L5/6", "L5/6-CC", "ext", "inb-pv", "inb-sst", "inb-sv2c", "inb-vip", "inb",
         "ast", "oli", "opc", "micro", "endo")
cs[, cell_type := factor(cell_type, levels = rev(ord))]
cs[, level := factor(level, levels = c("all genes", "length window"),
                     labels = c("All genes in set", "Shared length window only"))]
p2 <- ggplot(cs, aes(slope, cell_type, colour = set)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = slope - 1.96 * se, xmax = slope + 1.96 * se), height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  facet_wrap(~level) + scale_colour_manual(values = pal, labels = labs_set) +
  labs(x = "Slope of log2FC on log10(gene length)  (95% CI)", y = NULL,
       title = "Housekeeping genes show the length effect in every cell type",
       subtitle = "Neuron-specific genes are barely expressed in glia/endothelia, so their glial slopes are noisy")
save(p2, "fig2_length_slope_by_cell_type.png", 8, 5.5)

## Fig 3. sSNV rate by length decile ------------------------------------------------
dec <- fread("results/tables/mutation_rate_deciles.tsv")
r2  <- fread("results/tables/mutation_rate_decile_R2.tsv")
dec[, set_lab := labs_set[set]]
dec[, source := factor(source, levels = c("published (S16, 100 neurons)", "recount (13 donors, S13)"))]
p3 <- ggplot(dec, aes(decile, rel_rate, colour = set)) +
  geom_hline(yintercept = 1, colour = ink2, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.6) +
  geom_line(linewidth = 0.6, alpha = 0.6) + geom_point(size = 2.2) +
  facet_grid(source ~ set_lab) + scale_x_continuous(breaks = 1:10) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(x = "Gene-length decile (1 = shortest)", y = "sSNV rate per bp relative to set mean (exact Poisson 95% CI)",
       title = "The housekeeping length trend rests on the shortest decile",
       subtitle = sprintf("Housekeeping R^2 vs decile: published %.2f (p=%.3f) -> %.2f (p=%.2f) without decile 1",
                          r2[grepl("published", source) & set == "housekeeping", R2_all],
                          r2[grepl("published", source) & set == "housekeeping", p_all],
                          r2[grepl("published", source) & set == "housekeeping", R2_no_d1],
                          r2[grepl("published", source) & set == "housekeeping", p_no_d1]))
save(p3, "fig3_ssnv_rate_by_decile.png", 8, 6)

## Fig 4. RNA integrity vs age -------------------------------------------------------
d <- fread("results/tables/donor_covariates.tsv", colClasses = list(character = "case_id"))[!is.na(rin)]
cov <- fread("results/tables/donor_covariates_by_group.tsv")[variable == "rin"]
gpal <- c(infant = "#1baf7a", adult = "#2a78d6", elderly = "#eb6834")
p4 <- ggplot(d, aes(age, rin, colour = group)) +
  geom_smooth(data = d[group != "infant"], aes(group = 1), method = "lm", formula = y ~ x,
              colour = ink2, fill = grid, linewidth = 0.6) +
  geom_point(size = 3) + scale_colour_manual(values = gpal) +
  labs(x = "Donor age (years)", y = "RIN", title = "RNA integrity declines with donor age",
       subtitle = sprintf("Adult + elderly donors: Spearman rho = %.2f, p = %.3f (n = %d)",
                          cov$spearman_rho_with_age, cov$spearman_p, cov$n_adult + cov$n_elderly))
save(p4, "fig4_rin_vs_age.png", 6, 4.2)
message("figures written to results/figures/")
