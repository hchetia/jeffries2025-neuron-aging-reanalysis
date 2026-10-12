# One-time setup. Run in R/RStudio from the project root.
cran <- c("readxl", "data.table", "ggplot2", "patchwork", "MASS")
bioc <- c("org.Hs.eg.db", "GenomicRanges")
new <- setdiff(cran, rownames(installed.packages()))
if (length(new)) install.packages(new)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
newb <- setdiff(bioc, rownames(installed.packages()))
if (length(newb)) BiocManager::install(newb, update = FALSE)
