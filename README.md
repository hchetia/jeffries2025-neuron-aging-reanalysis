# Re-analysis of Jeffries et al. 2025 (*Nature*): gene length, ageing and somatic mutation in human neurons
Acknowledgement - Thanks to the original authors for depositing detailed supplementary data.
An independent re-test of the central claim of
[Jeffries et al., *Nature* 646:657 (2025)](https://doi.org/10.1038/s41586-025-09435-8):
that **short, highly expressed housekeeping genes accumulate somatic SNVs and lose expression
with age in human prefrontal-cortex neurons, while long neuron-specific genes are preserved.**

The study's raw data are controlled-access (dbGaP phs003445.v1.p1), so this re-analysis uses only
public material: the paper's supplementary tables, Ensembl/UCSC gene annotation, and two external
housekeeping-gene lists. Everything is in R and runs end-to-end in a couple of minutes.

## Summary

| Claim in the paper | Re-analysis | Verdict |
|---|---|---|
| Expression change in ageing rises with gene length in housekeeping genes (Fig. 5f; R² = 0.50) | Reproduced: R² = 0.51 (decile), 0.55 (continuous length), p < 10⁻²⁷⁰ | ✅ Reproduces |
| Housekeeping set is not an artefact of being defined on the same data | Same slope in genes also on an external HK list (0.56) and genes not on it (0.53) | ✅ Robust |
| Neuron-specific genes show a much weaker length relationship (Fig. 5g; R² = 0.20, P = 1.2 × 10⁻³) | R² = 0.24 but p = 6 × 10⁻¹⁹ at gene level; slope is 70% of HK in a shared length window |  Length effect is general, not housekeeping-specific (length and expression level not separable here) |
| sSNV rate per bp falls with gene length in housekeeping genes (Fig. 5i; R² = 0.44, P = 0.035) | Reproduced on deciles (R² = 0.44), but p = 0.12–0.20 without the shortest decile; HK vs neuron-specific slopes do not differ (p = 0.30–0.44); unchanged after removing likely-artifact recurrent calls | Fragile |
| Elderly and adult donors are comparable | RIN falls with age (ρ = −0.57, p = 0.021); depth, complexity, PMI, storage time balanced | Possible confound |

## Findings

### 1. The expression–length relationship is real, but it tracks gene length more than gene function
Using the per-gene fold changes in Table S11, the housekeeping result reproduces almost exactly
(excitatory neurons: R² = 0.51 with length decile as predictor vs 0.50 published).

The two gene classes barely overlap in length: neuron-specific genes are mostly long, and
housekeeping genes mostly short. Pooling both classes in excitatory neurons:

| Model for log2(elderly/adult) | R² |
|---|---|
| gene length only | **0.49** |
| gene class (HK vs neuron-specific) only | 0.02 |
| length + class | 0.50 |
| length × class | 0.53 |

Within the 9–108 kb window that both classes share, neuron-specific genes still have a clearly
positive slope (0.42 vs 0.59 for housekeeping; difference p = 0.008). The pattern holds in every
cell type, glia and endothelium included, though it is weaker there (Fig. 2). Removing ribosomal
and OXPHOS genes changes nothing (slope 0.55 vs 0.56).

**Interpretation:** short genes go down and long genes go up, in most cell types and across
gene classes. That fits the paper's mechanism, but it also fits any global length-dependent
effect, such as RNA integrity or intronic-read fraction with intron-inclusive counting.
The public tables can't separate these explanations (see Next steps).

**Caveat: length here also stands in for expression level.** Short housekeeping genes are also
the most highly expressed, and the public tables have no per-gene expression values to put in
the model. "Length explains 49%" therefore means length, or the expression level that goes with
it. The paper's own multiple regression (Table S15) included expression and found length had the
larger coefficient, but that can't be checked from the public tables.

![Expression vs length](results/figures/fig1_expression_vs_length.png)
![Length slope by cell type](results/figures/fig2_length_slope_by_cell_type.png)

### 2. The mutation–length result rests on one decile
The deposited per-SNV table (S13 "Detailed information") holds **19,585 of the 34,089** sSNVs
called across 100 neurons. Six of 19 donors are missing and four are partial, so Table S16
can't be rebuilt exactly. Recounting the available calls by coordinate overlap with GRCh37 gene
bodies gives decile counts about 0.6× the published ones, but the same shape.

| Housekeeping: rate vs length decile | R² | p |
|---|---|---|
| Published (S16, 100 neurons) | 0.44 | 0.035 |
| Published, decile 1 removed | 0.31 | 0.12 |
| Recount (13 donors) | 0.44 | 0.038 |
| Recount, decile 1 removed | 0.23 | 0.20 |

Gene-level negative-binomial models (`n_snv ~ log10(length) + offset(log length)`):

| Gene set | Rate change per 10× length | p |
|---|---|---|
| Housekeeping | 0.50× | 0.026 |
| Housekeeping, decile 1 removed | 0.64× | 0.19 |
| Housekeeping, 20–200 kb only | 1.26× | 0.64 |
| Neuron-specific | 0.73× | 0.28 |
| HK vs neuron-specific slope difference | — | **0.44** |
| All other protein-coding genes | 1.40× | 1 × 10⁻⁹ |

The elevated rate is real in the shortest ~180 housekeeping genes: in the published counts,
decile 1 has a 2.9× excess (exact Poisson 95% CI 1.9–4.2). It does not extend into a trend
across lengths, and it doesn't separate housekeeping from neuron-specific genes. Callable-region
size and GC content, which affect PTA/SCAN2 sensitivity and are not in the public data, could
also matter most for short GC-rich genes.

**Sensitivity to suspect calls.** Nine sites are called somatic in two or more unrelated donors.
That is very unlikely for true somatic mutations, so they are probably technical artifacts.
Another 468 rows repeat a site within one donor (clonal mutations shared by several neurons).
Neither changes the conclusions (`results/tables/mutation_sensitivity_recurrent_calls.tsv`):

| SNV calls used | Housekeeping p | Housekeeping p, decile 1 removed | HK vs NS slope p |
|---|---|---|---|
| All 19,585 | 0.026 | 0.19 | 0.44 |
| Cross-donor recurrent sites removed | 0.018 | 0.14 | 0.36 |
| + clonal duplicates collapsed | 0.019 | 0.13 | 0.30 |

![sSNV rate by decile](results/figures/fig3_ssnv_rate_by_decile.png)

### 3. RNA integrity declines with age in this cohort
Across the 16 adult and elderly donors with RIN, RIN correlates with age (Spearman ρ = −0.57,
p = 0.021; median 8.2 adult vs 6.7 elderly). Sequencing depth, genes per nucleus, PMI and storage
time do not differ. Library batches mix ages. The fold changes come from a cell-level Wilcoxon
test with no covariates, so RIN is not adjusted for in the published age effects.

![RIN vs age](results/figures/fig4_rin_vs_age.png)

## Limitations
- Per-gene fold changes are public only for the 1,802 housekeeping and 288 neuron-specific genes,
  so a genome-wide length test isn't possible from these tables.
- Gene length and expression level are confounded (short housekeeping genes are highly
  expressed), and per-gene expression isn't public, so the length effect can't be separated from
  an expression-level effect.
- The fold changes are the paper's own (Seurat, cell-level Wilcoxon). Without donor-level counts
  they can't be redone as pseudobulk or adjusted for RIN.
- The mutation analysis uses the 57% of sSNVs in the deposited table. Callable regions per neuron
  are not available, so rates are per annotated bp. A few calls look like artifacts (recurrent
  across unrelated donors); results are unchanged without them.

## Next steps
Repeat the expression test in an independent public cohort with donor-level counts and RIN/PMI,
such as Ling et al. 2024 (control BA46, ages 22–97) via NeMO. That means donor-level pseudobulk
(edgeR/limma) with RIN as a covariate, followed by a genome-wide test of fold change against
gene length *and* baseline expression, which would separate the two.

## Reproduce
```bash
Rscript scripts/install_packages.R   # one-time: readxl, data.table, ggplot2, patchwork, MASS, org.Hs.eg.db, GenomicRanges
bash run_all.sh                      # downloads public data, then runs scripts 02-08
```
Tested with R 4.3.3.

| Script | Purpose |
|---|---|
| `00_download.sh` | Supplementary files from Nature; GRCh37 gene table (annotables); external HK lists (scMerge) |
| `01_inspect_supp.R` | Inventory every sheet of the supplementary workbooks |
| `02_tidy_supp.R` | Workbooks → tidy TSVs |
| `03_gene_lengths.R` | Gene-body lengths (org.Hs.eg.db, hg38) and alias table |
| `04_reproduce_s16.R` | S13 completeness; gene-level sSNV counts on GRCh37; comparison with S16 |
| `05_mutation_vs_length.R` | Negative-binomial models and decile R² for sSNV rate vs length |
| `06_expression_vs_length.R` | Expression change vs length: reproduction, length-matched, by cell type, external HK lists |
| `07_donor_covariates.R` | RIN, PMI, storage, library QC vs age group |
| `08_figures.R` | Figures |

Outputs are in `results/tables/` and `results/figures/`. Genes whose annotated length disagrees
with the paper's decile by more than 5× (UCSC/Ensembl vs GENCODE gene-model differences) are excluded:
32 of 2,090 in the expression analysis and 6 in the mutation analysis.

## Data sources
- Jeffries, A.M. et al. Single-cell transcriptomic and genomic changes in the ageing human brain. *Nature* **646**, 657–666 (2025). Supplementary Tables 1–16.
- Ensembl GRCh37 gene coordinates via [`annotables`](https://github.com/stephenturner/annotables).
- Eisenberg, E. & Levanon, E.Y. Human housekeeping genes, revisited. *Trends Genet.* 29, 569 (2013); Lin, Y. et al. scMerge. *PNAS* 116, 9775 (2019). Lists via [`scMerge`](https://github.com/SydneyBioX/scMerge).
- `org.Hs.eg.db` (Bioconductor) for hg38 gene positions and symbol aliases.

*This is an independent re-analysis for learning purposes and is not affiliated with the original authors.*
