# Q4a. Rebuild the gene-level sSNV counts behind Table S16 from the per-SNV calls
# in Table S13, and compare with the published decile counts.
#
# SNV coordinates in S13 are GRCh37 (Annovar hg19), so gene bodies come from the
# Ensembl GRCh37 table shipped with the `annotables` R package.
# Run from project root:  Rscript scripts/04_reproduce_s16.R
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

e <- new.env(); load("data/raw/annotables/grch37.rda", envir = e)
g37 <- as.data.table(e$grch37)[chr %in% c(1:22, "X", "Y") & !is.na(symbol) & symbol != ""]
g37 <- g37[, .SD[which.max(end - start)], by = symbol]          # one locus per symbol
ali <- fread("data/processed/gene_aliases.tsv")
snv <- fread("data/processed/ssnv_calls.tsv", colClasses = list(character = "sample"))
summ <- fread("data/processed/neuron_summary.tsv", colClasses = list(character = "donor"))
gs  <- unique(fread("data/processed/geneset_log2fc.tsv")[, .(gene, length_decile, set)])

## 1. Completeness of the deposited per-SNV table ------------------------------
comp <- merge(summ[, .(neurons = .N, called_in_summary = sum(n_ssnv_called), age = age[1]), by = donor],
              snv[, .(rows_in_detailed = .N), by = .(donor = sample)], by = "donor", all = TRUE)
comp[is.na(rows_in_detailed), rows_in_detailed := 0]
setorder(comp, age)
fwrite(comp, "results/tables/s13_completeness_by_donor.tsv", sep = "\t")
message(sprintf("S13 detailed sheet: %d of %d called sSNVs (%.0f%%); %d of %d donors present",
                nrow(snv), sum(comp$called_in_summary), 100 * nrow(snv) / sum(comp$called_in_summary),
                sum(comp$rows_in_detailed > 0), nrow(comp)))

## 2. Gene-level counts by coordinate overlap -----------------------------------
gs[, symbol := fifelse(gene %in% g37$symbol, gene, ali$symbol[match(gene, ali$alias)])]
genes <- merge(g37[, .(symbol, chr, start, end, biotype)], gs[, .(symbol, gene, set, length_decile)],
               by = "symbol", all.x = TRUE)
genes <- genes[biotype == "protein_coding" | !is.na(set)]
genes[is.na(set), set := "other_protein_coding"]
genes[, gene_length := end - start + 1]

sgr <- GRanges(as.character(snv$Chr), IRanges(snv$Start, width = 1))
genes[, n_snv := countOverlaps(GRanges(chr, IRanges(start, end)), sgr)]
fwrite(genes, "data/processed/gene_level_snv_counts.tsv", sep = "\t")
message("gene-set genes with GRCh37 coordinates: ", genes[set != "other_protein_coding", .N], " / ", nrow(gs))

## 3. Compare with the published deciles ---------------------------------------
pub <- fread("data/processed/s16_decile_rates.tsv")
pub[, set := fifelse(grepl("House", set), "housekeeping", "neuron_specific")]
mine <- genes[set != "other_protein_coding", .(n_genes = .N, bp = sum(gene_length), mutations = sum(n_snv)),
              by = .(set, decile = length_decile)]
cmp <- merge(pub[, .(set, decile, pub_bp = total_bp, pub_mutations = mutations, pub_rate = rate)],
             mine, by = c("set", "decile"))
cmp[, ratio_pub_to_mine := round(pub_mutations / mutations, 2)]
setorder(cmp, set, decile)
fwrite(cmp, "results/tables/s16_reproduction.tsv", sep = "\t")
print(cmp)
