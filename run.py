"""Prepare results, list module commands, or explicitly execute one module."""
from pathlib import Path
from zipfile import ZipFile
import argparse
import csv
import io
import json
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def resolve(value):
    value = Path(value).expanduser()
    return value.resolve() if value.is_absolute() else (ROOT / value).resolve()

def prepare(config, profile):
    workspace = resolve(config['workspace'])
    if workspace.exists():
        raise SystemExit('Choose a new workspace; existing workspaces are preserved.')
    workspace.mkdir(parents=True)
    with ZipFile(ROOT / 'data/analysis_tables.zip') as archive:
        mapping = csv.DictReader(io.StringIO(archive.read('file_map.csv').decode('utf-8')))
        for row in mapping:
            name = row['workspace_file']
            if profile == 'analysis' and name not in (
                'inputs/bulk/samples.csv',
                'inputs/bulk/probe_mapping.csv',
                'enrichment/work/GENE_ID_MAPPING_DETAILS.csv'):
                continue
            destination = (workspace / name).resolve()
            if not destination.is_relative_to(workspace):
                raise SystemExit('Invalid relative input path.')
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(archive.read(row['archive_member']))
    shutil.copytree(ROOT / 'data/inputs/bulk', workspace / 'inputs/bulk', dirs_exist_ok=True)
    shutil.copyfile(ROOT / 'data/reference/bulk_gene_order.txt', workspace / 'inputs/bulk/genes.txt')
    for name in ('malignant','nonmalignant'):
        shutil.copytree(ROOT / 'data/inputs' / name, workspace / 'cellular' / name, dirs_exist_ok=True)
        original = workspace / 'cellular' / name / 'samples.csv'
        alias = 'SAMPLE_COLUMNS.csv' if name == 'malignant' else 'NONMALIGNANT_PSEUDOBULK_COLUMNS.csv'
        # Figure-data and library identifiers are shared across modules.
        shutil.copyfile(original, original.parent / alias)
    clinical = workspace / 'clinical'
    clinical.mkdir(exist_ok=True)
    for original, target in (
        ('data/reference/cohort_expression_manifest.json','clinical/SCORE_INPUT_CONFIG.json'),
        ('data/inputs/clinical/clinical_covariates.csv','clinical/CLINICAL_JOINED.csv'),
        ('data/reference/node_like_genes.csv','cohort_NODE_LIKE_SCORE_GENES.csv'),
        ('data/reference/clinical_sample_map.csv','clinical_TECHNICAL_SCORE_EXCLUSION_LEDGER.csv'),
        ('data/reference/development_folds.csv','models/inputs/DEVELOPMENT_FOLDS.csv')):
        destination = workspace / target
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT / original,destination)
    figures = resolve(config['figure_root'])
    figures.mkdir(parents=True,exist_ok=True)
    shutil.copytree(ROOT / 'data/figure_data',figures / 'source_data',dirs_exist_ok=True)
    for folder in ('main_figures','supplementary_figures'):
        (figures / folder).mkdir(exist_ok=True)
    print('Prepared '+profile+' workspace. No statistical analysis executed.')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','plan','run'))
    parser.add_argument('--profile',choices=('results','analysis'),default='results')
    parser.add_argument('--config',default='config/example.json')
    parser.add_argument('--stage')
    parser.add_argument('--execute',action='store_true')
    args = parser.parse_args()
    config = json.loads(resolve(args.config).read_text(encoding='utf-8'))
    if args.command == 'prepare':
        prepare(config,args.profile)
        return
    if args.stage not in config['stages']:
        raise SystemExit('Choose --stage from: '+', '.join(config['stages']))
    workspace = resolve(config['workspace'])
    commands = []
    for name in config['stages'][args.stage]:
        command = [config['python'] if name.endswith('.py') else config['rscript']]
        if name.endswith('.R'):
            command.append('--vanilla')
        command.append(str(ROOT / name))
        commands.append(command)
    print(json.dumps({'stage':args.stage,'commands':commands},indent=2))
    if args.command != 'run' or not args.execute:
        return
    if not workspace.is_dir():
        raise SystemExit('Prepare a workspace and restore the documented inputs first.')
    env = dict(os.environ,BCPN_REPOSITORY=str(ROOT),BCPN_WORKSPACE=str(workspace),
               BCPN_FIGURE_ROOT=str(resolve(config['figure_root'])))
    env['BCPN_R_LIBRARY'] = str(resolve(config['r_library'])) if config.get('r_library') else ''
    for command in commands:
        subprocess.run(command,cwd=ROOT,env=env,check=True)

if __name__ == '__main__':
    main()
