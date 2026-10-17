# Data guide

`tables/` is the main entry point for readers. The six supplementary workbooks correspond to the manuscript's Tables S1–S6 and include contents and notes. `Supplementary_Tables_CSV.zip` supplies full-precision text versions.

`figure_data/` contains the individual plotting inputs. `reference/` contains sample maps, selected-gene mappings, program definitions, reference scaling parameters, model coefficients and the fixed ridge reference. These files support computational reuse rather than adding new result tables to the manuscript.

`inputs/` contains the included processed bulk expression, selected patient-level pseudobulk counts and de-identified clinical covariates. File layouts and external resource requirements are described in `docs/INPUT_DATA.md`.

`analysis_tables.zip` supplies machine-readable intermediate tables for the analysis modules. Its `file_map.csv` uses relative workspace paths. The launcher can restore either precomputed results or analysis inputs into a separate workspace.
