from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "manual_curation_outputs_merged" / "figure_data" / "appendix_tables"

MODEL_STAGE_PATHS: List[Tuple[str, Path]] = [
    ("iDT1294Photo", ROOT / "iDT1294Photo.json"),
    ("DSM123", ROOT / "Models" / "DSM123.json"),
    ("updated_consensus", ROOT / "updated_consensus.json"),
    ("purple_bacteriav_DSM123", ROOT / "Models" / "purple_bacteriav_DSM123.json"),
]

BBH_PARSED = ROOT / "bbh" / "Rpal_BisA53_vs_DSM123_parsed.csv"
REACTION_GENE_REL = ROOT / "reaction_gene_relationships.csv"
MANUAL_REPLACE = ROOT / "Models" / "old" / "DSM123_blast_gene_rename_output_rpal_final" / "applied_manual_replacements.csv"
MANUAL_DELETED_GENES = (
    ROOT / "Models" / "old" / "DSM123_blast_gene_rename_output_rpal_final" / "applied_manual_deleted_genes.csv"
)

PID_THRESHOLD = 65.0
COV_THRESHOLD = 0.2
EVALUE_THRESHOLD = 1e-3

LOGIC_WORDS = {"and", "or", "not"}
GENE_TOKEN_PATTERN = re.compile(r"[A-Za-z0-9_.:-]+")


def write_csv(path: Path, rows: List[Dict[str, object]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_model_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def extract_gene_ids(model_json: Dict[str, object]) -> Set[str]:
    out: Set[str] = set()
    for gene in model_json.get("genes", []):
        if not isinstance(gene, dict):
            continue
        gid = gene.get("id")
        if isinstance(gid, str) and gid:
            out.add(gid)
    return out


def extract_reaction_ids(model_json: Dict[str, object]) -> Set[str]:
    out: Set[str] = set()
    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        rid = rxn.get("id")
        if isinstance(rid, str) and rid:
            out.add(rid)
    return out


def extract_metabolite_ids(model_json: Dict[str, object]) -> Set[str]:
    out: Set[str] = set()
    for met in model_json.get("metabolites", []):
        if not isinstance(met, dict):
            continue
        mid = met.get("id")
        if isinstance(mid, str) and mid:
            out.add(mid)
    return out


def gene_prefix(gene_id: str) -> str:
    if "_" in gene_id:
        return gene_id.split("_", 1)[0]
    return gene_id


def format_counter(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return " | ".join(f"{k}:{v}" for k, v in counter.most_common())


def parse_genes_from_gpr(gpr: str) -> Set[str]:
    if not gpr:
        return set()
    out: Set[str] = set()
    for token in GENE_TOKEN_PATTERN.findall(gpr):
        if token.lower() in LOGIC_WORDS:
            continue
        out.add(token)
    return out


def build_gene_to_reactions(model_json: Dict[str, object]) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = defaultdict(set)
    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        rid = rxn.get("id")
        if not isinstance(rid, str) or not rid:
            continue
        genes = parse_genes_from_gpr(str(rxn.get("gene_reaction_rule", "") or ""))
        for gid in genes:
            mapping[gid].add(rid)
    return mapping


def sorted_join(items: Iterable[str], sep: str = " | ") -> str:
    values = sorted({x for x in items if x})
    return sep.join(values)


def transition_summary_rows(
    stage_models: Dict[str, Dict[str, object]],
) -> List[Dict[str, object]]:
    stage_names = [name for name, _ in MODEL_STAGE_PATHS]
    rows: List[Dict[str, object]] = []
    for i in range(len(stage_names) - 1):
        old_name = stage_names[i]
        new_name = stage_names[i + 1]
        old_model = stage_models[old_name]
        new_model = stage_models[new_name]

        old_genes = extract_gene_ids(old_model)
        new_genes = extract_gene_ids(new_model)
        old_rxns = extract_reaction_ids(old_model)
        new_rxns = extract_reaction_ids(new_model)
        old_mets = extract_metabolite_ids(old_model)
        new_mets = extract_metabolite_ids(new_model)

        genes_added = sorted(new_genes - old_genes)
        genes_removed = sorted(old_genes - new_genes)
        rxns_added = sorted(new_rxns - old_rxns)
        rxns_removed = sorted(old_rxns - new_rxns)
        mets_added = sorted(new_mets - old_mets)
        mets_removed = sorted(old_mets - new_mets)

        rows.append(
            {
                "transition": f"{old_name} -> {new_name}",
                "old_model": old_name,
                "new_model": new_name,
                "old_model_path": str(dict(MODEL_STAGE_PATHS)[old_name]),
                "new_model_path": str(dict(MODEL_STAGE_PATHS)[new_name]),
                "genes_old": len(old_genes),
                "genes_new": len(new_genes),
                "genes_added": len(genes_added),
                "genes_removed": len(genes_removed),
                "genes_common": len(old_genes & new_genes),
                "gene_prefix_added_breakdown": format_counter(Counter(gene_prefix(g) for g in genes_added)),
                "gene_prefix_removed_breakdown": format_counter(Counter(gene_prefix(g) for g in genes_removed)),
                "reactions_old": len(old_rxns),
                "reactions_new": len(new_rxns),
                "reactions_added": len(rxns_added),
                "reactions_removed": len(rxns_removed),
                "reactions_common": len(old_rxns & new_rxns),
                "metabolites_old": len(old_mets),
                "metabolites_new": len(new_mets),
                "metabolites_added_by_id": len(mets_added),
                "metabolites_removed_by_id": len(mets_removed),
                "metabolites_common_by_id": len(old_mets & new_mets),
            }
        )
    return rows


def transition_compact_list_rows(stage_models: Dict[str, Dict[str, object]]) -> List[Dict[str, object]]:
    stage_names = [name for name, _ in MODEL_STAGE_PATHS]
    rows: List[Dict[str, object]] = []
    for i in range(len(stage_names) - 1):
        old_name = stage_names[i]
        new_name = stage_names[i + 1]
        old_model = stage_models[old_name]
        new_model = stage_models[new_name]

        old_genes = extract_gene_ids(old_model)
        new_genes = extract_gene_ids(new_model)
        old_rxns = extract_reaction_ids(old_model)
        new_rxns = extract_reaction_ids(new_model)
        old_mets = extract_metabolite_ids(old_model)
        new_mets = extract_metabolite_ids(new_model)

        rows.append(
            {
                "transition": f"{old_name} -> {new_name}",
                "added_gene_ids": sorted_join(new_genes - old_genes),
                "removed_gene_ids": sorted_join(old_genes - new_genes),
                "added_reaction_ids": sorted_join(new_rxns - old_rxns),
                "removed_reaction_ids": sorted_join(old_rxns - new_rxns),
                "added_metabolite_ids_by_id": sorted_join(new_mets - old_mets),
                "removed_metabolite_ids_by_id": sorted_join(old_mets - new_mets),
            }
        )
    return rows


def load_best_bbh_tables(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not path.exists():
        empty = pd.DataFrame(columns=["gene", "subject", "PID", "eVal", "COV", "BBH"])
        return empty, empty

    df = pd.read_csv(path)
    if df.empty:
        empty = pd.DataFrame(columns=["gene", "subject", "PID", "eVal", "COV", "BBH"])
        return empty, empty

    order = df.sort_values(
        by=["PID", "eVal", "COV", "bitScore"],
        ascending=[False, True, False, False],
        kind="stable",
    )
    best_by_gene = order.groupby("gene", as_index=False).first()
    best_by_subject = order.groupby("subject", as_index=False).first()
    return best_by_gene, best_by_subject


def load_reaction_gene_relationships(path: Path) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = defaultdict(set)
    if not path.exists():
        return mapping
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        rid = str(row.get("Reaction", "")).strip()
        genes_raw = str(row.get("Genes", "")).strip()
        if not rid or not genes_raw or genes_raw.lower() == "nan":
            continue
        for token in genes_raw.split(","):
            gid = token.strip()
            if gid:
                mapping[gid].add(rid)
    return mapping


def load_manual_replacements(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not path.exists():
        return mapping
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        old_gene = str(row.get("old_gene_id", "")).strip()
        new_gene = str(row.get("new_gene_id", "")).strip()
        if old_gene and new_gene and old_gene.lower() != "nan" and new_gene.lower() != "nan":
            mapping[old_gene] = new_gene
    return mapping


def load_manual_deleted_genes(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not path.exists():
        return mapping
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        gid = str(row.get("gene_id", "")).strip()
        reason = str(row.get("delete_reason", "")).strip()
        if gid and gid.lower() != "nan":
            mapping[gid] = reason
    return mapping


def reaction_info_for_gene(
    gene_id: str,
    old_gene_to_rxns: Dict[str, Set[str]],
    new_gene_to_rxns: Dict[str, Set[str]],
) -> Tuple[int, int, str, str]:
    old_rxns = old_gene_to_rxns.get(gene_id, set())
    new_rxns = new_gene_to_rxns.get(gene_id, set())
    return len(old_rxns), len(new_rxns), sorted_join(old_rxns), sorted_join(new_rxns)


def build_transition_gene_change_rows(
    stage_models: Dict[str, Dict[str, object]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []

    # Transition 1: iDT1294Photo -> DSM123
    t1_old_name = "iDT1294Photo"
    t1_new_name = "DSM123"
    old_genes_t1 = extract_gene_ids(stage_models[t1_old_name])
    new_genes_t1 = extract_gene_ids(stage_models[t1_new_name])
    t1_old_rxn_map = build_gene_to_reactions(stage_models[t1_old_name])
    t1_new_rxn_map = build_gene_to_reactions(stage_models[t1_new_name])

    best_by_gene, best_by_subject = load_best_bbh_tables(BBH_PARSED)
    best_by_gene_map = {str(row["gene"]): row for _, row in best_by_gene.iterrows()}
    best_by_subject_map = {str(row["subject"]): row for _, row in best_by_subject.iterrows()}

    replacement_targets_t1: Set[str] = set()
    for old_gene in sorted(old_genes_t1):
        hit = best_by_gene_map.get(old_gene)
        if hit is None:
            old_count, new_count, old_rxns, new_rxns = reaction_info_for_gene(old_gene, t1_old_rxn_map, t1_new_rxn_map)
            rows.append(
                {
                    "transition": f"{t1_old_name} -> {t1_new_name}",
                    "change_mode": "remove",
                    "old_gene_id": old_gene,
                    "new_gene_id": "",
                    "reason_category": "no_bbh_hit_remove_reference_gene",
                    "reason_detail": "No BBH hit found in parsed BLAST table; reference gene removed in DSM123 draft.",
                    "reason_source": str(BBH_PARSED),
                    "confidence": "direct",
                    "bbh_pid": "",
                    "bbh_evalue": "",
                    "bbh_cov": "",
                    "old_reaction_count": old_count,
                    "new_reaction_count": new_count,
                    "old_reaction_ids": old_rxns,
                    "new_reaction_ids": new_rxns,
                }
            )
            continue

        subject = str(hit.get("subject", "") or "").strip()
        pid = float(hit.get("PID")) if pd.notna(hit.get("PID")) else None
        evalue = float(hit.get("eVal")) if pd.notna(hit.get("eVal")) else None
        cov = float(hit.get("COV")) if pd.notna(hit.get("COV")) else None
        passed = (
            pid is not None
            and evalue is not None
            and cov is not None
            and (pid >= PID_THRESHOLD)
            and (cov >= COV_THRESHOLD)
            and (evalue <= EVALUE_THRESHOLD)
        )

        old_count, _, old_rxns, _ = reaction_info_for_gene(old_gene, t1_old_rxn_map, t1_new_rxn_map)

        if subject and subject in new_genes_t1:
            replacement_targets_t1.add(subject)
            _, new_count, _, new_rxns = reaction_info_for_gene(subject, t1_old_rxn_map, t1_new_rxn_map)
            category = "bbh_ortholog_replacement_pass_threshold" if passed else "bbh_ortholog_replacement_selected"
            detail = (
                "Reference gene replaced by DSM123 ortholog using BBH mapping and threshold filter."
                if passed
                else "Reference gene replaced by DSM123 ortholog using BBH mapping; retained in draft model."
            )
            rows.append(
                {
                    "transition": f"{t1_old_name} -> {t1_new_name}",
                    "change_mode": "replace",
                    "old_gene_id": old_gene,
                    "new_gene_id": subject,
                    "reason_category": category,
                    "reason_detail": detail,
                    "reason_source": f"{BBH_PARSED}; thresholds(pid>={PID_THRESHOLD}, cov>={COV_THRESHOLD}, e<={EVALUE_THRESHOLD})",
                    "confidence": "direct",
                    "bbh_pid": "" if pid is None else pid,
                    "bbh_evalue": "" if evalue is None else evalue,
                    "bbh_cov": "" if cov is None else cov,
                    "old_reaction_count": old_count,
                    "new_reaction_count": new_count,
                    "old_reaction_ids": old_rxns,
                    "new_reaction_ids": new_rxns,
                }
            )
        else:
            rows.append(
                {
                    "transition": f"{t1_old_name} -> {t1_new_name}",
                    "change_mode": "remove",
                    "old_gene_id": old_gene,
                    "new_gene_id": "",
                    "reason_category": "bbh_hit_not_selected_gene_removed",
                    "reason_detail": "BBH hit exists but mapped target gene not retained in DSM123 draft after model pruning.",
                    "reason_source": str(BBH_PARSED),
                    "confidence": "inferred",
                    "bbh_pid": "" if pid is None else pid,
                    "bbh_evalue": "" if evalue is None else evalue,
                    "bbh_cov": "" if cov is None else cov,
                    "old_reaction_count": old_count,
                    "new_reaction_count": 0,
                    "old_reaction_ids": old_rxns,
                    "new_reaction_ids": "",
                }
            )

    for new_gene in sorted(new_genes_t1 - replacement_targets_t1):
        hit = best_by_subject_map.get(new_gene)
        old_candidate = str(hit.get("gene", "")).strip() if hit is not None else ""
        pid = float(hit.get("PID")) if (hit is not None and pd.notna(hit.get("PID"))) else None
        evalue = float(hit.get("eVal")) if (hit is not None and pd.notna(hit.get("eVal"))) else None
        cov = float(hit.get("COV")) if (hit is not None and pd.notna(hit.get("COV"))) else None
        _, new_count, _, new_rxns = reaction_info_for_gene(new_gene, t1_old_rxn_map, t1_new_rxn_map)
        rows.append(
            {
                "transition": f"{t1_old_name} -> {t1_new_name}",
                "change_mode": "add",
                "old_gene_id": old_candidate if old_candidate in old_genes_t1 else "",
                "new_gene_id": new_gene,
                "reason_category": "target_gene_present_without_unique_replacement_pair",
                "reason_detail": "DSM123 gene retained in draft but not consumed by one-to-one replacement pairing from reference genes.",
                "reason_source": str(BBH_PARSED),
                "confidence": "inferred",
                "bbh_pid": "" if pid is None else pid,
                "bbh_evalue": "" if evalue is None else evalue,
                "bbh_cov": "" if cov is None else cov,
                "old_reaction_count": 0,
                "new_reaction_count": new_count,
                "old_reaction_ids": "",
                "new_reaction_ids": new_rxns,
            }
        )

    # Transition 2: DSM123 -> updated_consensus
    t2_old_name = "DSM123"
    t2_new_name = "updated_consensus"
    old_genes_t2 = extract_gene_ids(stage_models[t2_old_name])
    new_genes_t2 = extract_gene_ids(stage_models[t2_new_name])
    t2_old_rxn_map = build_gene_to_reactions(stage_models[t2_old_name])
    t2_new_rxn_map = build_gene_to_reactions(stage_models[t2_new_name])
    added_t2 = sorted(new_genes_t2 - old_genes_t2)
    rel_map = load_reaction_gene_relationships(REACTION_GENE_REL)

    for new_gene in added_t2:
        added_rxns = rel_map.get(new_gene, set())
        _, new_count, _, new_rxns = reaction_info_for_gene(new_gene, t2_old_rxn_map, t2_new_rxn_map)
        if added_rxns:
            category = "flux_supported_gapfill_gene_added"
            detail = "Gene introduced with flux-supported missing reaction import from reference model."
            source = f"{REACTION_GENE_REL}; missing_reactions.csv"
            confidence = "direct"
        else:
            category = "gapfill_related_gene_added_unmatched_table"
            detail = "Gene appears in updated model but was not found in reaction_gene_relationships table."
            source = str(REACTION_GENE_REL)
            confidence = "inferred"

        rows.append(
            {
                "transition": f"{t2_old_name} -> {t2_new_name}",
                "change_mode": "add",
                "old_gene_id": "",
                "new_gene_id": new_gene,
                "reason_category": category,
                "reason_detail": detail,
                "reason_source": source,
                "confidence": confidence,
                "bbh_pid": "",
                "bbh_evalue": "",
                "bbh_cov": "",
                "old_reaction_count": 0,
                "new_reaction_count": new_count,
                "old_reaction_ids": "",
                "new_reaction_ids": sorted_join(new_rxns.split(" | ") if new_rxns else added_rxns),
            }
        )

    # Transition 3: updated_consensus -> purple_bacteriav_DSM123
    t3_old_name = "updated_consensus"
    t3_new_name = "purple_bacteriav_DSM123"
    old_genes_t3 = extract_gene_ids(stage_models[t3_old_name])
    new_genes_t3 = extract_gene_ids(stage_models[t3_new_name])
    t3_old_rxn_map = build_gene_to_reactions(stage_models[t3_old_name])
    t3_new_rxn_map = build_gene_to_reactions(stage_models[t3_new_name])
    added_t3 = set(new_genes_t3 - old_genes_t3)
    removed_t3 = set(old_genes_t3 - new_genes_t3)

    manual_replace_map = load_manual_replacements(MANUAL_REPLACE)
    manual_deleted_map = load_manual_deleted_genes(MANUAL_DELETED_GENES)

    consumed_added: Set[str] = set()
    consumed_removed: Set[str] = set()

    # Direct manual replacement pairs
    for old_gene, new_gene in sorted(manual_replace_map.items()):
        if old_gene not in old_genes_t3 or new_gene not in new_genes_t3:
            continue
        old_count = len(t3_old_rxn_map.get(old_gene, set()))
        new_count = len(t3_new_rxn_map.get(new_gene, set()))
        old_rxns = sorted_join(t3_old_rxn_map.get(old_gene, set()))
        new_rxns = sorted_join(t3_new_rxn_map.get(new_gene, set()))

        rows.append(
            {
                "transition": f"{t3_old_name} -> {t3_new_name}",
                "change_mode": "replace",
                "old_gene_id": old_gene,
                "new_gene_id": new_gene,
                "reason_category": "manual_replacement_from_curator_table",
                "reason_detail": "Old gene ID replaced by curated new gene ID from manual replacement table.",
                "reason_source": str(MANUAL_REPLACE),
                "confidence": "direct",
                "bbh_pid": "",
                "bbh_evalue": "",
                "bbh_cov": "",
                "old_reaction_count": old_count,
                "new_reaction_count": new_count,
                "old_reaction_ids": old_rxns,
                "new_reaction_ids": new_rxns,
            }
        )

        if old_gene in removed_t3:
            consumed_removed.add(old_gene)
        if new_gene in added_t3:
            consumed_added.add(new_gene)

    # Explicit manual gene deletions should be recorded before inferred pairing.
    for old_gene in sorted(manual_deleted_map.keys()):
        if old_gene not in removed_t3 or old_gene in consumed_removed:
            continue
        old_set = t3_old_rxn_map.get(old_gene, set())
        rows.append(
            {
                "transition": f"{t3_old_name} -> {t3_new_name}",
                "change_mode": "remove",
                "old_gene_id": old_gene,
                "new_gene_id": "",
                "reason_category": "explicit_manual_gene_deletion",
                "reason_detail": f"Gene explicitly deleted by curator instruction: {manual_deleted_map[old_gene]}.",
                "reason_source": str(MANUAL_DELETED_GENES),
                "confidence": "direct",
                "bbh_pid": "",
                "bbh_evalue": "",
                "bbh_cov": "",
                "old_reaction_count": len(old_set),
                "new_reaction_count": 0,
                "old_reaction_ids": sorted_join(old_set),
                "new_reaction_ids": "",
            }
        )
        consumed_removed.add(old_gene)

    # Inferred replacement by exact reaction-set transfer
    remaining_removed = sorted(removed_t3 - consumed_removed)
    remaining_added = sorted(added_t3 - consumed_added)

    added_by_rxnset: Dict[Tuple[str, ...], List[str]] = defaultdict(list)
    for new_gene in remaining_added:
        rxn_set = tuple(sorted(t3_new_rxn_map.get(new_gene, set())))
        if rxn_set:
            added_by_rxnset[rxn_set].append(new_gene)

    for old_gene in remaining_removed:
        old_set = tuple(sorted(t3_old_rxn_map.get(old_gene, set())))
        if not old_set:
            continue
        candidates = added_by_rxnset.get(old_set, [])
        if len(candidates) != 1:
            continue
        new_gene = candidates[0]
        if new_gene in consumed_added:
            continue
        rows.append(
            {
                "transition": f"{t3_old_name} -> {t3_new_name}",
                "change_mode": "replace",
                "old_gene_id": old_gene,
                "new_gene_id": new_gene,
                "reason_category": "inferred_replacement_identical_reaction_set",
                "reason_detail": "Old/new gene carry identical reaction assignments across transition, indicating GPR relabeling.",
                "reason_source": "set_difference + GPR reaction-set identity",
                "confidence": "inferred",
                "bbh_pid": "",
                "bbh_evalue": "",
                "bbh_cov": "",
                "old_reaction_count": len(old_set),
                "new_reaction_count": len(old_set),
                "old_reaction_ids": sorted_join(old_set),
                "new_reaction_ids": sorted_join(old_set),
            }
        )
        consumed_removed.add(old_gene)
        consumed_added.add(new_gene)

    # Remaining removed genes
    for old_gene in sorted(removed_t3 - consumed_removed):
        old_set = t3_old_rxn_map.get(old_gene, set())
        if old_gene in manual_deleted_map:
            category = "explicit_manual_gene_deletion"
            detail = f"Gene explicitly deleted by curator instruction: {manual_deleted_map[old_gene]}."
            source = str(MANUAL_DELETED_GENES)
            confidence = "direct"
        elif old_gene.startswith("RPE_"):
            category = "reference_gene_cleanup_in_final_model"
            detail = "Reference-prefixed gene removed during final DSM123-focused cleanup."
            source = "updated vs final gene set difference"
            confidence = "inferred"
        else:
            category = "curation_pruning_gene_removed"
            detail = "Gene removed during final manual curation and GPR cleanup."
            source = "updated vs final gene set difference"
            confidence = "inferred"

        rows.append(
            {
                "transition": f"{t3_old_name} -> {t3_new_name}",
                "change_mode": "remove",
                "old_gene_id": old_gene,
                "new_gene_id": "",
                "reason_category": category,
                "reason_detail": detail,
                "reason_source": source,
                "confidence": confidence,
                "bbh_pid": "",
                "bbh_evalue": "",
                "bbh_cov": "",
                "old_reaction_count": len(old_set),
                "new_reaction_count": 0,
                "old_reaction_ids": sorted_join(old_set),
                "new_reaction_ids": "",
            }
        )

    # Remaining added genes
    for new_gene in sorted(added_t3 - consumed_added):
        new_set = t3_new_rxn_map.get(new_gene, set())
        prefix = gene_prefix(new_gene)
        if prefix in {"PGIDNB", "b1210"}:
            category = "non_target_prefix_gene_added_or_retained"
            detail = "Gene with non-ACXYSJ prefix appears in final model GPR after manual curation."
            source = "updated vs final gene set difference"
            confidence = "inferred"
        else:
            category = "target_gene_added_in_final_curation"
            detail = "DSM123-prefixed gene added in final model during GPR/manual refinement."
            source = "updated vs final gene set difference"
            confidence = "inferred"

        rows.append(
            {
                "transition": f"{t3_old_name} -> {t3_new_name}",
                "change_mode": "add",
                "old_gene_id": "",
                "new_gene_id": new_gene,
                "reason_category": category,
                "reason_detail": detail,
                "reason_source": source,
                "confidence": confidence,
                "bbh_pid": "",
                "bbh_evalue": "",
                "bbh_cov": "",
                "old_reaction_count": 0,
                "new_reaction_count": len(new_set),
                "old_reaction_ids": "",
                "new_reaction_ids": sorted_join(new_set),
            }
        )

    rows.sort(
        key=lambda x: (
            str(x["transition"]),
            str(x["change_mode"]),
            str(x["reason_category"]),
            str(x["old_gene_id"]),
            str(x["new_gene_id"]),
        )
    )
    return rows


def summarize_reason_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counts: Counter[Tuple[str, str, str]] = Counter()
    for row in rows:
        key = (str(row["transition"]), str(row["change_mode"]), str(row["reason_category"]))
        counts[key] += 1
    out: List[Dict[str, object]] = []
    for (transition, change_mode, reason_category), count in sorted(counts.items()):
        out.append(
            {
                "transition": transition,
                "change_mode": change_mode,
                "reason_category": reason_category,
                "row_count": count,
            }
        )
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    stage_models = {name: load_model_json(path) for name, path in MODEL_STAGE_PATHS}

    detail_rows = build_transition_gene_change_rows(stage_models)
    reasons_path = OUT_DIR / "appendix_table_transition_gene_change_reasons.csv"
    write_csv(
        reasons_path,
        detail_rows,
        [
            "transition",
            "change_mode",
            "old_gene_id",
            "new_gene_id",
            "reason_category",
            "reason_detail",
            "reason_source",
            "confidence",
            "bbh_pid",
            "bbh_evalue",
            "bbh_cov",
            "old_reaction_count",
            "new_reaction_count",
            "old_reaction_ids",
            "new_reaction_ids",
        ],
    )

    for path in OUT_DIR.glob("appendix_table_*"):
        if path.name == reasons_path.name:
            continue
        if path.is_file():
            try:
                path.unlink()
            except PermissionError:
                pass

    print("Exported appendix table:")
    print(f"- transition_gene_change_reasons: {reasons_path}")


if __name__ == "__main__":
    main()
