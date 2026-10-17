"""Export the author-prepared graphical abstract and optionally render its PDF."""

from pathlib import Path
import argparse
import os
import shutil


REPOSITORY = Path(__file__).resolve().parents[2]
SOURCE = REPOSITORY / 'figures'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default=os.environ.get('BCPN_FIGURE_ROOT', 'workspace/figures'))
    parser.add_argument('--render-raster', action='store_true',
                        help='Render new PNG/TIFF exports from the supplied PDF instead of copying them.')
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error('--dpi must be positive')
    output = Path(args.output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    for extension in ['svg', 'pdf', 'png', 'tiff']:
        source = SOURCE / f'GraphicalAbstract.{extension}'
        target = output / f'GRAPHICAL_ABSTRACT.{extension}'
        if not source.is_file():
            raise SystemExit(f'Required deposited artwork is missing: {source.name}')
        if source.resolve() != target:
            shutil.copyfile(source, target)
    source = SOURCE / 'GraphicalAbstract_caption.txt'
    target = output / 'GRAPHICAL_ABSTRACT_CAPTION.txt'
    if source.resolve() != target:
        shutil.copyfile(source, target)
    if args.render_raster:
        import pypdfium2 as pdfium
        with pdfium.PdfDocument(str(output / 'GRAPHICAL_ABSTRACT.pdf')) as document:
            page = document[0]
            try:
                picture = page.render(scale=args.dpi/72).to_pil().convert('RGB')
                picture.save(output / 'GRAPHICAL_ABSTRACT.png', dpi=(args.dpi,args.dpi))
                picture.save(output / 'GRAPHICAL_ABSTRACT.tiff', compression='tiff_lzw', dpi=(args.dpi,args.dpi))
            finally:
                page.close()
    print(f'Graphical abstract exported to {output}.')


if __name__ == '__main__':
    main()
