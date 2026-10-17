"""Reuse the selected METABRIC binary matrix and identity ledger; no outcomes."""
from pathlib import Path
import os
import json
import numpy as np
import pandas as pd

root = Path(os.environ['BCPN_WORKSPACE']).resolve()
out = root / 'clinical_context'
assert not (out / 'work/METABRIC_PREPARATION.json').exists()
config = next(c for c in json.loads((root / 'clinical/SCORE_INPUT_CONFIG.json').read_text()) if c['cohort'] == 'METABRIC')
assert config['patients'] == 1964 and config['genes'] == 19349
binary = root / config['binary_file']
assert binary.stat().st_size == config['patients'] * config['genes'] * 8
matrix = np.memmap(binary, dtype='<f8', mode='r', shape=(config['patients'], config['genes']))
defs = pd.read_csv(out / 'work/NODAL_SIGNATURE_DEFINITIONS.csv')
required = sorted(set(g for members in defs.full_members for g in members.split(';')))
indices = {g: i for i, g in enumerate(config['background_genes'])}
available = [g for g in required if g in indices]
values = np.asarray(matrix[:, [indices[g] for g in available]])
assert np.isfinite(values).all()
means, sds = values.mean(axis=0), values.std(axis=0, ddof=1)
usable = np.isfinite(sds) & (sds > 0)
genes = np.asarray(available)[usable].tolist()
values = values[:, usable]
means, sds = means[usable], sds[usable]
z = (values - means) / sds
pd.DataFrame(z.T, index=genes, columns=config['patient_ids']).rename_axis('gene').to_csv(out / 'work/METABRIC_SIGNATURE_GENE_Z.csv')
pd.DataFrame(dict(gene=genes,cohort_gene_mean=means,cohort_gene_SD=sds)).to_csv(out / 'work/METABRIC_GENE_STANDARDIZATION.csv', index=False)
coverage = []
for r in defs.itertuples():
    members = r.full_members.split(';')
    used = sorted(set(members) & set(genes))
    coverage.append(dict(score=r.score, full_set_N=len(members), used_N=len(used), coverage_pct=100*len(used)/len(members),
                         missing_or_constant=';'.join(sorted(set(members)-set(used))), used_members=';'.join(used),
                         normalization='Across all 1964 previously technically eligible primary-tumor samples; sample SD; no nodal-label conditioning'))
pd.DataFrame(coverage).to_csv(out / 'work/METABRIC_SIGNATURE_COVERAGE.csv', index=False)
clinical = pd.read_csv(Path(os.environ['BCPN_REPOSITORY']) / 'data/reference/metabric_mask.csv',dtype=str)
assert clinical.SAMPLE_ID.tolist()==config['patient_ids']
assert len(clinical)==1964 and clinical.PATIENT_ID.is_unique and clinical.SAMPLE_TYPE.eq('Primary').all()
count = pd.to_numeric(clinical.LYMPH_NODES_EXAMINED_POSITIVE, errors='coerce')
valid = count.notna() & (count >= 0) & np.isclose(count.fillna(0), count.fillna(0).round())
clinical['positive_node_count'] = count
clinical['nodal_count_eligible'] = valid
clinical['nodal_group'] = np.select([valid & count.eq(0),valid & count.between(1,3),valid & count.ge(4)],['N0','LOW','HIGH'],default='NOT_EVALUABLE')
clinical['exclusion_reason'] = np.where(valid,'','MISSING_OR_INVALID_POSITIVE_NODE_COUNT')
clinical.to_csv(out / 'work/METABRIC_CLINICAL_MAP.csv', index=False)
