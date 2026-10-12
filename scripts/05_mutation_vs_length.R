# Q4b. Does the per-bp sSNV rate fall with gene length in housekeeping genes but
# not in neuron-specific genes (paper Fig. 5i,j)? Tested gene-by-gene with
# negative-binomial models instead of a 10-point regression on deciles.
#
#   n_snv ~ NB(mu),  log(mu) = a + b * log10(length) + log(length)
#   b < 0  => rate per bp decreases with gene length (rate ~ length^b)
# Run from project root:  Rscript scripts/05_mutation_vs_length.R
suppressPackageStartupMessages({ library(data.table); library(MASS) })
g   <- fread("data/processed/gene_level_snv_counts.tsv")
pub <- fread("data/processed/s16_decile_rates.tsv")
g[, l10 := log10(gene_length)]
# Same consistency filter as script 06: drop paper-set genes whose GRCh37 length is >5x off their decile
g[set != "other_protein_coding", dec_med := as.numeric(median(gene_length)), by = .(set, length_decile)]
message("gene-set genes dropped for length/decile disagreement: ",
        g[set != "other_protein_coding" & (gene_length > 5 * dec_med | gene_length < dec_med / 5), .N])
g <- g[set == "other_protein_coding" | (gene_length <= 5 * dec_med & gene_length >= dec_med / 5)]

fit_b <- function(d, label) {
  m <- tryCatch(glm.nb(n_snv ~ l10 + offset(log(gene_length)), data = d, control = glm.control(maxit = 100)),
                error = function(e) glm(n_snv ~ l10 + offset(log(gene_length)), family = quasipoisson, data = d))
  co <- summary(m)$coefficients["l10", ]
  data.table(model = label, n_genes = nrow(d), n_snv = sum(d$n_snv),
             b = co[1], se = co[2], p = co[4], family = class(m)[1])
}

res <- rbind(
  fit_b(g[set == "housekeeping"],        "housekeeping (paper set)"),
  fit_b(g[set == "neuron_specific"],     "neuron-specific (paper set)"),
  fit_b(g[set == "other_protein_coding"], "all other protein-coding genes"),
  fit_b(g[set == "housekeeping" & length_decile > 1], "housekeeping, decile 1 removed")
)

# Same-length window: is HK different from NS once length ranges overlap?
win <- g[set %in% c("housekeeping", "neuron_specific") & gene_length >= 2e4 & gene_length <= 2e5]
res <- rbind(res,
  fit_b(win[set == "housekeeping"],    "housekeeping, 20-200 kb only"),
  fit_b(win[set == "neuron_specific"], "neuron-specific, 20-200 kb only"))
mi <- glm.nb(n_snv ~ set * l10 + offset(log(gene_length)), data = g[set %in% c("housekeeping", "neuron_specific")])
inter <- summary(mi)$coefficients["setneuron_specific:l10", ]
res[, rate_ratio_per_10x_length := round(10^b, 2)]
print(res, digits = 3)
message(sprintf("HK vs NS difference in slope (interaction): %.3f (SE %.3f), p = %.3g", inter[1], inter[2], inter[4]))
fwrite(res, "results/tables/mutation_rate_vs_length_models.tsv", sep = "\t")

## Decile-level: published values and our recount, with exact Poisson CIs ------
dec <- g[set != "other_protein_coding", .(mutations = sum(n_snv), bp = as.numeric(sum(gene_length))), by = .(set, decile = length_decile)]
dec[, rel_rate := (mutations / bp) / (sum(mutations) / sum(bp)), by = set]
dec[, c("lo", "hi") := as.data.table(t(mapply(function(x, e) poisson.test(x, e)$conf.int, mutations,
                                                bp * sum(mutations) / sum(bp)))), by = set]
dec[, source := "recount (13 donors, S13)"]
pd <- pub[, .(set = fifelse(grepl("House", set), "housekeeping", "neuron_specific"), decile, mutations, bp = as.numeric(total_bp))]
pd[, rel_rate := (mutations / bp) / (sum(mutations) / sum(bp)), by = set]
pd[, c("lo", "hi") := as.data.table(t(mapply(function(x, e) poisson.test(x, e)$conf.int, mutations,
                                               bp * sum(mutations) / sum(bp)))), by = set]
pd[, source := "published (S16, 100 neurons)"]
dec <- rbind(dec, pd)
fwrite(dec, "results/tables/mutation_rate_deciles.tsv", sep = "\t")

# The paper's statistic: Pearson R^2 of rate vs decile, with and without decile 1
r2 <- dec[, .(R2_all = cor(rel_rate, decile)^2, p_all = cor.test(rel_rate, decile)$p.value,
              R2_no_d1 = cor(rel_rate[decile > 1], decile[decile > 1])^2,
              p_no_d1 = cor.test(rel_rate[decile > 1], decile[decile > 1])$p.value), by = .(source, set)]
print(r2, digits = 3)
fwrite(r2, "results/tables/mutation_rate_decile_R2.tsv", sep = "\t")

## Sensitivity: recurrent and clonal calls -------------------------------------
# A site called somatic in two unrelated donors is very unlikely to be a true somatic
# mutation (probable artifact); a site repeated within one donor is a clonal mutation
# shared by several neurons. Refit the key models with these removed / collapsed.
suppressPackageStartupMessages(library(GenomicRanges))
snv <- fread("data/processed/ssnv_calls.tsv", colClasses = list(character = "sample"))
cross <- snv[, uniqueN(sample), by = .(Chr, Start)][V1 > 1]
calls <- list(
  "all calls"                         = snv,
  "cross-donor recurrent sites removed" = snv[!cross, on = .(Chr, Start)],
  "+ clonal duplicates collapsed"     = unique(snv[!cross, on = .(Chr, Start)], by = c("Chr", "Start", "sample")))
ggr <- GRanges(g$chr, IRanges(g$start, g$end))
sens <- rbindlist(lapply(names(calls), function(nm) {
  s <- calls[[nm]]
  g[, n_snv := countOverlaps(ggr, GRanges(as.character(s$Chr), IRanges(s$Start, width = 1)))]
  hk  <- fit_b(g[set == "housekeeping"], "HK")
  hk1 <- fit_b(g[set == "housekeeping" & length_decile > 1], "HK no d1")
  mi  <- glm.nb(n_snv ~ set * l10 + offset(log(gene_length)), data = g[set %in% c("housekeeping", "neuron_specific")])
  data.table(calls = nm, n_calls = nrow(s), hk_b = hk$b, hk_p = hk$p, hk_no_d1_b = hk1$b, hk_no_d1_p = hk1$p,
             hk_vs_ns_p = summary(mi)$coefficients["setneuron_specific:l10", 4])
}))
message("cross-donor recurrent sites: ", nrow(cross), "; within-donor duplicate rows: ",
        nrow(snv[!cross, on = .(Chr, Start)]) - nrow(calls[[3]]))
print(sens, digits = 3)
fwrite(sens, "results/tables/mutation_sensitivity_recurrent_calls.tsv", sep = "\t")
