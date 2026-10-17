# Analysis methods

## Paired bulk discovery

GSE57968 contributes 36 patients with one primary tumor and one involved regional node per patient. Deposited processed log2 expression is analyzed using the supplied unambiguous selected-probe map and 12,502-gene universe. Probe choice used the highest all-specimen mean, without tissue labels. The design is `~ patient + tissue`, with Primary as the reference and the `tissueNode` coefficient estimating Node minus Primary. limma empirical Bayes moderation uses `robust=FALSE, trend=FALSE`. All tested genes share one Benjamini–Hochberg family; candidates have q < 0.05 without a fold-change threshold.

## Enrichment and tissue scores

Directional upper-tail hypergeometric tests use database-annotated tested genes as backgrounds and include eligible zero-overlap terms. Correction is separate by database and direction. Ranked enrichment uses the complete moderated-t ranking; Hallmark, GO, Reactome and KEGG retain separate correction families. The Hallmark call uses seed 20261005, minSize 15, maxSize 500, eps 0, sampleSize 101, nPermSimple 1000, scoreType std, gseaParam 1 and nproc 1. Other ranked collections use seed 20261006 and eps 1e-10 with the remaining parameters unchanged.

ESTIMATE and MCP-counter readouts use two-sided paired Wilcoxon tests with `exact=FALSE`, continuity correction and zero differences omitted from signed ranks. Correction families contain four and ten slots, respectively. CD8 T-cell scoring is non-evaluable where its required marker is unavailable. Scores are inferred expression readouts, not measured fractions.

## Cellular expression support

Patients are the biological replicates. The exploratory 20-cell rule retains five paired patients; the initial 50-cell rule retained two. Per-compartment TMM normalization uses the full 25,301-feature background. Filtering uses `filterByExpr` with full raw library sizes, min.count 10, min.total.count 15, large.n 10 and min.prop 0.7. Negative-binomial quasi-likelihood nuisance fitting uses the complete filtered background, with candidate-only tests and the original fitted quantities retained. `glmQLFit` uses legacy FALSE, abundance.trend TRUE, robust FALSE and winsor.tail.p c(0.05,0.10); `glmQLFTest` uses tissueNode and poisson.bound FALSE.

Malignant tests share a candidate-only correction family. All evaluable nonmalignant candidate-by-population tests share one joint correction family. Concordant q < 0.05 supports an expression change within a population; it does not quantify that population's contribution to the bulk change. Non-evaluable and unsupported candidates are retained. Descriptive patient-resampling intervals use all 3,125 ordered five-patient samples and type-8 percentiles. Threshold and leave-one-patient-out summaries do not refit inferential models.

## Molecular readouts

See `MOLECULAR_READOUTS.md` for formulas, missing features, pooled transformations and primary-defined reference strata. The PAM50 available-feature implementation was adopted after continuous-readout analysis but before examining subtype assignments. The 23 continuous readouts share one correction family, and the three relative research scores form another. Correlated readouts are not independent clinical assay validations.

## Clinical context and supplementary models

Five METABRIC program scores use gene means and sample standard deviations from all 1,964 expression-eligible primary tumors. The 1,890 patients with valid positive-node counts enter the unadjusted association analyses. Spearman correlations, Kruskal–Wallis comparisons and exploratory ordinal associations retain separate five-test families. Proportional odds was not independently assessed.

Supplementary global-score and prediction results retain their original clinical coding, reference parameters and comparison families. The compressed-model winner used METABRIC performance; associated model comparisons are exploratory and selection-conditioned. These analyses did not establish incremental clinical usefulness. Full estimates, including nonsignificant results, are supplied in Table S6.
