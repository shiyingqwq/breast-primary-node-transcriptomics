# Input data

## Included inputs

- `data/inputs/bulk/GSE57968_expression.txt.gz`: the GPL571 expression table extracted from the deposited processed matrix. GEO contact metadata are omitted; expression values and sample identifiers are preserved. Values are RMA/log2; no additional normalization is required.
- `data/reference/bulk_samples.csv` and `bulk_pairs.csv`: explicit GSM, patient and tissue identities for 36 pairs.
- `data/reference/bulk_gene_order.txt`: order of the 12,502 tested genes. Selected probes and annotation mappings are supplied in the machine-readable tables.
- `data/inputs/malignant/` and `nonmalignant/`: GSE167036 patient-level raw-count pseudobulks, gene names and sample maps. `counts.float64` uses little-endian 64-bit values in row-major order, with genes as rows and sample-map rows as columns. Counts are integer-valued. The malignant matrix has 25,301 rows and 10 columns; the nonmalignant matrix has 25,301 rows and 90 columns across nine eligible populations. Mast cells remain represented as non-evaluable in the result tables.
- `data/inputs/clinical/clinical_covariates.csv`: de-identified analysis covariates and source endpoint coding. These are supplementary-analysis inputs, not new prospective clinical records.
- `data/reference/cohort_expression_manifest.json`: ordered clinical-cohort sample and gene identifiers, dimensions and relative expression-file locations. The large clinical expression matrices are not included.

## External resources

1. Obtain MSigDB 2025.1.Hs Hallmark, GO biological-process and Reactome GMT files from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb). The workspace paths are `resources/h.all.v2025.1.Hs.symbols.gmt`, `enrichment/sources/GO_BP.v2025.1.Hs.gmt` and `enrichment/sources/REACTOME.v2025.1.Hs.gmt`. Preserve the supplied resource terms.
2. KEGG raw mappings are not redistributed. The analyses used the human pathway snapshot retrieved on 6 October 2026, with a pathway release date of 5 October 2026. Restore `enrichment/sources/KEGG_GENE_PATHWAY.tsv` and `KEGG_PATHWAYS.tsv` under the applicable [KEGG terms](https://www.kegg.jp/kegg/legal.html). A current query need not reproduce that snapshot. Published enrichment outputs are supplied for inspection.
3. Install genefu 2.40.0 for its released PAM50 robust model, classification/rescaling functions and signature definition files. The code uses the documented 43 available PAM50 features. No third-party package binaries are included.
4. Restore the MCP-counter author implementation and marker files as `enrichment/sources/MCPcounter.R`, `MCP_genes.tsv` and `MCP_probesets.tsv`. ESTIMATE scoring uses tidyestimate 1.1.1.
5. Restore processed METABRIC and GSE96058 expression after the documented median duplicate-gene aggregation. Each binary is a little-endian float64 matrix in sample-major order, using exactly the sample and 19,349-gene orders in `cohort_expression_manifest.json`. Locations are `clinical/expression/METABRIC.float64` and `clinical/expression/GSE96058.float64` within the chosen workspace. Original source links are in the main README.
6. The supplementary compressed-model modules additionally require sample-major expression z-score CSVs in `models/inputs/`, with first column `patient_id` and the 869 foreground genes in the supplied reference-model feature order: `GSE96058_EXPRESSION_Z.csv`, `METABRIC_EXPRESSION_Z.csv` and `GSE96058_CLINICAL_TRAINING.csv`. Cohort-wide gene means and sample standard deviations must be used before selecting nodal-evaluable patients. Development folds and the fixed ridge reference are supplied. These optional modules must be explicitly requested through the launcher.

## Workspace preparation

`run.py prepare --profile results` unpacks all machine-readable tables using `data/analysis_tables.zip/file_map.csv`. `--profile analysis` retains only sample, probe and annotation inputs, and copies the included processed input files. Reference resources are restored to relative workspace paths. No workstation-specific directories or private access credentials are required by the distribution.
