# Gene-body length (min start .. max end over all transcripts) for every gene symbol,
# from org.Hs.eg.db (UCSC hg38 gene positions). Writes data/processed/gene_lengths.tsv
# Run from the project root:  Rscript scripts/03_gene_lengths.R
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(data.table) })

s <- as.data.table(toTable(org.Hs.egCHRLOC));    setnames(s, c("entrez", "start", "chr"))
e <- as.data.table(toTable(org.Hs.egCHRLOCEND)); setnames(e, c("entrez", "end", "chr"))
loc <- cbind(s, end = e$end)                      # rows are paired transcript start/end
loc <- loc[chr %in% c(1:22, "X", "Y")]
loc[, strand := fifelse(start < 0, "-", "+")]
loc[, `:=`(start = abs(start), end = abs(end))]
# Some genes have placements on >1 chromosome (PAR, duplicated loci): keep the
# chromosome with the most transcript records before taking the span.
loc[, n_chr := .N, by = .(entrez, chr)]
loc <- loc[, .SD[n_chr == max(n_chr) & chr == chr[which.max(n_chr)]], by = entrez]
len <- loc[, .(chr = chr[1], start = min(start), end = max(end), strand = strand[1]), by = entrez]
len[, gene_length := end - start + 1]

sym <- as.data.table(toTable(org.Hs.egSYMBOL)); setnames(sym, c("entrez", "symbol"))
len <- merge(len, sym, by = "entrez")

# Alias table so older symbols in the paper's tables still match
ali <- as.data.table(toTable(org.Hs.egALIAS2EG)); setnames(ali, c("entrez", "alias"))
ali <- ali[!alias %in% sym$symbol]                  # never let an alias shadow an official symbol
ali <- ali[, if (.N == 1) .SD, by = alias]          # drop ambiguous aliases
ali <- merge(ali, len[, .(entrez, symbol)], by = "entrez")

fwrite(len[, .(symbol, entrez, chr, start, end, strand, gene_length)], "data/processed/gene_lengths.tsv", sep = "\t")
fwrite(ali[, .(alias, symbol)], "data/processed/gene_aliases.tsv", sep = "\t")
message("wrote gene_lengths.tsv (", nrow(len), " genes) and gene_aliases.tsv (", nrow(ali), ")")
