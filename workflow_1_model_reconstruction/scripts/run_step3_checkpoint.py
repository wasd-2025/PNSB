from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import appdirs
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_cache_dir = ROOT / '.cobra_cache'
_cache_dir.mkdir(parents=True, exist_ok=True)
appdirs.user_cache_dir = lambda *args, **kwargs: str(_cache_dir.resolve())

import cobra

MODEL_PATH = ROOT / 'Models' / 'DSM123_manual_working.json'
OUTPUT_DIR = ROOT / 'manual_curation_outputs_merged' / 'step3_review'


def normalize_name(name: str) -> str:
    if not name:
        return ''
    text = name.lower()
    text = re.sub(r'[^a-z0-9]+', ' ', text)
    return re.sub(r'\s+', ' ', text).strip()


def build_df(rows: List[Dict[str, object]], columns: List[str], sort_by: List[str]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        return df
    return df.sort_values(by=sort_by, kind='stable')


def run_fba_checkpoint(model: cobra.Model) -> pd.DataFrame:
    solution = model.optimize()
    objective_value = 0.0
    if solution.objective_value is not None:
        objective_value = float(solution.objective_value)

    record = {
        'model_file': str(MODEL_PATH),
        'status': solution.status,
        'objective_value': objective_value,
        'objective_direction': model.objective.direction,
        'objective_expression': str(model.objective.expression),
        'reaction_count': len(model.reactions),
        'metabolite_count': len(model.metabolites),
        'gene_count': len(model.genes),
    }
    return pd.DataFrame([record])


def audit_metabolite_consistency(model: cobra.Model) -> Dict[str, pd.DataFrame]:
    suffix_rows: List[Dict[str, object]] = []
    mismatch_rows: List[Dict[str, object]] = []
    orphan_rows: List[Dict[str, object]] = []
    missing_formula_rows: List[Dict[str, object]] = []

    grouped: Dict[tuple[str, str, str], List[cobra.Metabolite]] = defaultdict(list)

    for met in model.metabolites:
        met_id = met.id
        met_name = met.name or ''
        met_formula = met.formula or ''
        met_comp = met.compartment or ''

        if met_id.endswith('_u'):
            suffix_rows.append(
                {
                    'metabolite_id': met_id,
                    'name': met_name,
                    'formula': met_formula,
                    'compartment': met_comp,
                    'reaction_count': len(met.reactions),
                    'reaction_ids': ' | '.join(sorted(r.id for r in met.reactions)),
                }
            )

        match = re.search(r'_([A-Za-z0-9]+)$', met_id)
        id_suffix = match.group(1) if match else ''
        if id_suffix and met_comp and id_suffix != met_comp:
            mismatch_rows.append(
                {
                    'metabolite_id': met_id,
                    'name': met_name,
                    'formula': met_formula,
                    'id_suffix': id_suffix,
                    'compartment': met_comp,
                    'reaction_count': len(met.reactions),
                }
            )

        if not met_formula:
            missing_formula_rows.append(
                {
                    'metabolite_id': met_id,
                    'name': met_name,
                    'compartment': met_comp,
                    'reaction_count': len(met.reactions),
                }
            )

        boundary_rxns = sorted(r.id for r in met.reactions if r.boundary)
        internal_rxns = sorted(r.id for r in met.reactions if not r.boundary)

        if len(met.reactions) == 0:
            orphan_type = 'no_reaction'
        elif len(internal_rxns) == 0:
            orphan_type = 'boundary_only'
        elif len(internal_rxns) == 1:
            orphan_type = 'single_internal_reaction'
        else:
            orphan_type = ''

        if orphan_type:
            orphan_rows.append(
                {
                    'orphan_type': orphan_type,
                    'metabolite_id': met_id,
                    'name': met_name,
                    'formula': met_formula,
                    'compartment': met_comp,
                    'total_reaction_count': len(met.reactions),
                    'internal_reaction_count': len(internal_rxns),
                    'boundary_reaction_count': len(boundary_rxns),
                    'internal_reactions': ' | '.join(internal_rxns),
                    'boundary_reactions': ' | '.join(boundary_rxns),
                }
            )

        if met_formula and met_comp:
            grouped[(met_comp, met_formula, normalize_name(met_name))].append(met)

    duplicate_rows: List[Dict[str, object]] = []
    for (compartment, formula, normalized_name), mets in grouped.items():
        if len(mets) <= 1:
            continue
        ids = sorted(m.id for m in mets)
        names = sorted({m.name or '' for m in mets})
        duplicate_rows.append(
            {
                'group_key': f'{compartment}|{formula}|{normalized_name}',
                'compartment': compartment,
                'formula': formula,
                'normalized_name': normalized_name,
                'member_count': len(mets),
                'metabolite_ids': ' | '.join(ids),
                'metabolite_names': ' | '.join(names),
            }
        )

    return {
        'suffix_u_remaining': build_df(
            suffix_rows,
            ['metabolite_id', 'name', 'formula', 'compartment', 'reaction_count', 'reaction_ids'],
            ['metabolite_id'],
        ),
        'compartment_mismatch': build_df(
            mismatch_rows,
            ['metabolite_id', 'name', 'formula', 'id_suffix', 'compartment', 'reaction_count'],
            ['metabolite_id'],
        ),
        'orphans': build_df(
            orphan_rows,
            [
                'orphan_type',
                'metabolite_id',
                'name',
                'formula',
                'compartment',
                'total_reaction_count',
                'internal_reaction_count',
                'boundary_reaction_count',
                'internal_reactions',
                'boundary_reactions',
            ],
            ['orphan_type', 'metabolite_id'],
        ),
        'missing_formula': build_df(
            missing_formula_rows,
            ['metabolite_id', 'name', 'compartment', 'reaction_count'],
            ['metabolite_id'],
        ),
        'duplicate_candidates': build_df(
            duplicate_rows,
            [
                'group_key',
                'compartment',
                'formula',
                'normalized_name',
                'member_count',
                'metabolite_ids',
                'metabolite_names',
            ],
            ['group_key'],
        ),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    model = cobra.io.load_json_model(str(MODEL_PATH))

    fba_df = run_fba_checkpoint(model)
    audits = audit_metabolite_consistency(model)

    fba_path = OUTPUT_DIR / 'step3_table_fba_checkpoint.csv'
    suffix_path = OUTPUT_DIR / 'step3_table_suffix_u_remaining.csv'
    mismatch_path = OUTPUT_DIR / 'step3_table_compartment_mismatch.csv'
    orphan_path = OUTPUT_DIR / 'step3_table_orphan_metabolites.csv'
    missing_formula_path = OUTPUT_DIR / 'step3_table_missing_formula.csv'
    duplicate_path = OUTPUT_DIR / 'step3_table_duplicate_metabolite_candidates.csv'

    fba_df.to_csv(fba_path, index=False, encoding='utf-8-sig')
    audits['suffix_u_remaining'].to_csv(suffix_path, index=False, encoding='utf-8-sig')
    audits['compartment_mismatch'].to_csv(mismatch_path, index=False, encoding='utf-8-sig')
    audits['orphans'].to_csv(orphan_path, index=False, encoding='utf-8-sig')
    audits['missing_formula'].to_csv(missing_formula_path, index=False, encoding='utf-8-sig')
    audits['duplicate_candidates'].to_csv(duplicate_path, index=False, encoding='utf-8-sig')

    summary = {
        'model': str(MODEL_PATH),
        'counts': {
            'fba_records': int(len(fba_df)),
            'suffix_u_remaining': int(len(audits['suffix_u_remaining'])),
            'compartment_mismatch': int(len(audits['compartment_mismatch'])),
            'orphan_metabolites': int(len(audits['orphans'])),
            'missing_formula': int(len(audits['missing_formula'])),
            'duplicate_groups': int(len(audits['duplicate_candidates'])),
        },
        'tables': {
            'fba_checkpoint': str(fba_path),
            'suffix_u_remaining': str(suffix_path),
            'compartment_mismatch': str(mismatch_path),
            'orphan_metabolites': str(orphan_path),
            'missing_formula': str(missing_formula_path),
            'duplicate_metabolite_candidates': str(duplicate_path),
        },
    }
    (OUTPUT_DIR / 'step3_review_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

    print('Generated step3 review artifacts.')
    print(f'step3: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
