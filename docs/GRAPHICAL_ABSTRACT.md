# Graphical abstract

The graphical abstract is supplied as author-prepared editable SVG artwork, a vector PDF, and PNG/TIFF exports. The SVG canvas has a 2:1 aspect ratio; the PNG and TIFF files contain 4252 by 2126 pixels at 600 dpi. Arial is used for the SVG text, and fonts are embedded in the PDF. Install Arial or review font substitution when editing the SVG on another system.

Anatomical, cellular and expression-profile drawings are schematic. The abstract summarizes previously reported results; its preparation and export do not compute new biological results.

## Caption

Graphical abstract. Patient-aware analysis of 36 matched primary breast tumor and regional lymph node pairs identified 890 differential genes and coordinated stromal, B-cell receptor and cell-cycle program differences. Research PAM50 assignments differed in 8/36 pairs (22.2%), including six luminal A-to-B transitions. Independent cellular profiles from five paired patients supported selected malignant-cell changes; 1,890 METABRIC primary tumors provided unadjusted nodal-burden context. Anatomical, cellular and expression-profile illustrations are schematic. Research molecular readouts are not calibrated clinical assays. Artwork prepared by the authors.

## Editing and exports

Open `figures/GraphicalAbstract.svg` in an SVG editor to edit the artwork. Review typography, clipping and image rights after making changes, and regenerate corresponding exports through the editor.

To copy the reported artwork files to the configured output directory:

```sh
python run.py run --stage graphical-abstract --execute
```

Alternatively, export to a chosen directory directly:

```sh
python src/figures/graphical_abstract.py --output-dir workspace/figures
```

The following command renders PNG and TIFF files from the supplied vector PDF at 600 dpi:

```sh
python src/figures/graphical_abstract.py --output-dir workspace/figures --render-raster --dpi 600
```

Raster dimensions may differ by a few pixels from the supplied SVG-derived exports because of PDF page rounding; antialiasing can also differ between rendering engines. This command copies the supplied SVG and PDF and does not regenerate the PDF from the SVG.
