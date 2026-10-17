# Reproduction instructions

## Environment

Bulk and cellular analyses used R 4.5.3 and Bioconductor 3.21 with limma 3.64.3, edgeR 4.6.3 and fgsea 1.34.2. Molecular classification used genefu 2.40.0. Additional package versions are listed in `environment/`. Preserve the relevant module's environment rather than assuming all analyses used one identical library.

Python figure modules require NumPy, pandas, matplotlib and Pillow. Optional raster rendering of the graphical abstract uses pypdfium2 and Pillow. Install the listed Python dependencies in your chosen environment:

```sh
python -m pip install -r requirements.txt
```

Edit `config/example.json` to select your executables, optional R library and workspace. Paths are relative to the repository unless an explicit user configuration supplies an absolute path.

## Replot reported results

```sh
python run.py prepare --profile results
python run.py plan --stage figures
python run.py run --stage figures --execute
```

Preparation copies precomputed tables and figure values into `workspace/`; it refuses to overwrite an existing workspace. Plotting computes display geometry, not new scores or statistical tests. The separately supplied reported figure PDFs remain the publication record.

The graphical abstract has an author-prepared editable SVG source and accompanying PDF, PNG and TIFF files in `figures/`. Its export module copies these supplied files into the selected figure output directory. Optional raster rendering from its PDF is described in [Graphical abstract](GRAPHICAL_ABSTRACT.md); it does not recreate the author's SVG layout or perform statistical analysis.

## Re-execute analyses

Use a new configuration and workspace. Prepare only inputs to avoid mixing saved results with newly generated outputs:

```sh
python run.py prepare --profile analysis --config config/example.json
python run.py plan --stage bulk --config config/example.json
python run.py run --stage bulk --config config/example.json --execute
```

The order for the paired analysis is `bulk`, `enrichment`, `molecular`, then `clinical-context`. Cellular modules are `malignant` and `nonmalignant`, with the bulk candidate table available first. Supplementary clinical scoring and prediction modules have additional inputs described in `INPUT_DATA.md`. Run only the explicitly chosen module; the launcher does not automatically train prediction models.

Restore version-matched gene-set collections, third-party implementations and processed clinical expression before dependent modules. Missing resources must be restored from their documented sources; current source downloads can differ from the versions used in the study. Do not substitute gene sets, probes, patients or reference scaling based on resulting performance.

This reorganized distribution has not undergone a complete end-to-end rerun. Results supplied here are the previously computed study results. The material supports inspection, figure reproduction and method reuse; exact analytical reproduction additionally requires the specified inputs and environment.
