from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import appdirs
import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
_cache_dir = ROOT / ".cobra_cache"
_cache_dir.mkdir(parents=True, exist_ok=True)
appdirs.user_cache_dir = lambda *args, **kwargs: str(_cache_dir.resolve())

import cobra

BASELINE_MODEL = ROOT / "Models" / "DSM123.json"
FINAL_MODEL = ROOT / "Models" / "DSM123_manual_working.json"
REFERENCE_MODEL = ROOT / "iDT1294Photo.json"
BBH_PARSED = ROOT / "bbh" / "Rpal_BisA53_vs_DSM123_parsed.csv"
TARGET_GB = ROOT / "genomes" / "DSM123.gb"

OUTPUT_ROOT = ROOT / "manual_curation_outputs_merged"
STEP1_DIR = OUTPUT_ROOT / "step1_review"
STEP2_DIR = OUTPUT_ROOT / "step2_review"

OLD_BRIDGE_TABLE = STEP2_DIR / "step2_table_compartment_bridge_reactions.csv"
BRIDGE_ID_PATTERNS = (
    re.compile(r"^TR_.+_(cp|ep|ce)(?:_\d+)?$", re.IGNORECASE),
    re.compile(r"^BRIDGE_(CP|EP|CE)_.+$", re.IGNORECASE),
)

IDENTITY_THRESHOLD = 50.0
EVALUE_THRESHOLD = 1e-3


def configure_cobra_cache() -> None:
    cache_dir = ROOT / ".cobra_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    appdirs.user_cache_dir = lambda *args, **kwargs: str(cache_dir.resolve())


def normalize_name(name: str) -> str:
    if not name:
        return ""
    text = name.lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def canonical_formula(formula: str | None) -> str:
    if not formula:
        return ""
    return re.sub(r"\s+", "", str(formula)).upper()


def get_metabolite_compartment(met: cobra.Metabolite) -> str:
    comp = (met.compartment or "").strip()
    if comp:
        return comp
    match = re.search(r"_([A-Za-z0-9]+)$", met.id)
    return match.group(1) if match else ""


def reaction_equation_from_json(
    reaction_json: Dict[str, object],
    metabolite_name_map: Dict[str, str],
) -> str:
    mets = reaction_json.get("metabolites", {})
    if not isinstance(mets, dict):
        return ""
    left: List[str] = []
    right: List[str] = []
    for met_id, coeff in mets.items():
        try:
            val = float(coeff)
        except (TypeError, ValueError):
            continue
        met_name = metabolite_name_map.get(met_id, met_id)
        token = f"{abs(val):g} {met_id} ({met_name})" if abs(val) != 1 else f"{met_id} ({met_name})"
        if val < 0:
            left.append(token)
        elif val > 0:
            right.append(token)
    arrow = "<=>" if reaction_json.get("lower_bound", 0) < 0 else "-->"
    return f"{' + '.join(left)} {arrow} {' + '.join(right)}".strip()


def parse_gene_annotations(genbank_path: Path) -> Dict[str, Dict[str, str]]:
    annotations: Dict[str, Dict[str, str]] = {}
    with genbank_path.open("r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                locus_tag = ""
                if "locus_tag" in feature.qualifiers and feature.qualifiers["locus_tag"]:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers and feature.qualifiers["gene"]:
                    locus_tag = feature.qualifiers["gene"][0]
                if not locus_tag:
                    continue
                annotations[locus_tag] = {
                    "symbol": feature.qualifiers.get("gene", [""])[0] if feature.qualifiers.get("gene") else "",
                    "product": feature.qualifiers.get("product", [""])[0] if feature.qualifiers.get("product") else "",
                }
    return annotations


def make_reaction_index(model_json: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    for rxn in model_json.get("reactions", []):
        rid = rxn.get("id")
        if isinstance(rid, str):
            out[rid] = rxn
    return out


def make_metabolite_name_map(model_json: Dict[str, object]) -> Dict[str, str]:
    return {
        m.get("id", ""): m.get("name", "")
        for m in model_json.get("metabolites", [])
        if isinstance(m.get("id"), str)
    }


def list_gene_reactions(model: cobra.Model) -> Dict[str, List[str]]:
    mapping: Dict[str, List[str]] = {}
    for gene in model.genes:
        mapping[gene.id] = sorted(r.id for r in gene.reactions)
    return mapping


def load_bridge_ids_from_table() -> set[str]:
    if not OLD_BRIDGE_TABLE.exists():
        return set()
    try:
        df = pd.read_csv(OLD_BRIDGE_TABLE)
    except Exception:
        return set()
    if "reaction_id" not in df.columns:
        return set()
    return {str(x) for x in df["reaction_id"].dropna().tolist()}


def is_bridge_reaction_id(reaction_id: str) -> bool:
    return any(pattern.match(reaction_id) for pattern in BRIDGE_ID_PATTERNS)


def collect_bridge_reaction_ids(
    current_final: cobra.Model,
    baseline_rxn_ids: set[str],
) -> set[str]:
    bridge_ids = load_bridge_ids_from_table()
    for rxn in current_final.reactions:
        if rxn.id in baseline_rxn_ids:
            continue
        if is_bridge_reaction_id(rxn.id):
            bridge_ids.add(rxn.id)
            continue
        if (rxn.name or "").lower().startswith("auto bridge "):
            bridge_ids.add(rxn.id)
    return bridge_ids


def add_metabolite_from_template(
    model: cobra.Model,
    template: cobra.Metabolite,
) -> cobra.Metabolite:
    if template.id in model.metabolites:
        return model.metabolites.get_by_id(template.id)
    new_met = cobra.Metabolite(
        template.id,
        name=template.name,
        formula=template.formula,
        charge=template.charge,
        compartment=template.compartment,
    )
    model.add_metabolites([new_met])
    return new_met


def reconstruct_pre_step2_model() -> cobra.Model:
    baseline = cobra.io.load_json_model(str(BASELINE_MODEL))
    current_final = cobra.io.load_json_model(str(FINAL_MODEL))
    baseline_rxn_map = {rxn.id: rxn for rxn in baseline.reactions}
    bridge_ids = collect_bridge_reaction_ids(current_final, set(baseline_rxn_map))

    if bridge_ids:
        to_remove = [rxn for rxn in current_final.reactions if rxn.id in bridge_ids]
        if to_remove:
            current_final.remove_reactions(to_remove, remove_orphans=False)

    for rxn in list(current_final.reactions):
        if rxn.id not in baseline_rxn_map:
            continue
        baseline_rxn = baseline_rxn_map[rxn.id]
        rxn.subtract_metabolites(dict(rxn.metabolites))
        stoich = {}
        for met, coeff in baseline_rxn.metabolites.items():
            restored = add_metabolite_from_template(current_final, met)
            stoich[restored] = coeff
        rxn.add_metabolites(stoich)

    unused = [m for m in current_final.metabolites if len(m.reactions) == 0]
    if unused:
        current_final.remove_metabolites(unused, destructive=False)
    return current_final


def run_u_suffix_cleanup_with_audit(model: cobra.Model) -> Tuple[pd.DataFrame, Dict[str, int]]:
    records: List[Dict[str, object]] = []
    renamed = 0
    merged = 0

    for met in list(model.metabolites):
        old_id = met.id
        if not old_id.endswith("_u"):
            continue
        new_id = old_id[:-2]
        old_comp = met.compartment
        old_name = met.name
        old_formula = met.formula
        affected_rxns = sorted(r.id for r in met.reactions)

        if new_id in model.metabolites:
            target = model.metabolites.get_by_id(new_id)
            for rxn in list(met.reactions):
                coeff = rxn.metabolites[met]
                rxn.add_metabolites({target: coeff, met: -coeff}, combine=True)
            model.remove_metabolites([met], destructive=False)
            merged += 1
            action = "merged_u_into_existing"
            new_comp = target.compartment
        else:
            met.id = new_id
            comp_match = re.search(r"_([a-zA-Z0-9]+)$", new_id)
            if comp_match:
                met.compartment = comp_match.group(1)
            renamed += 1
            action = "renamed_u"
            new_comp = met.compartment

        records.append(
            {
                "action": action,
                "old_metabolite_id": old_id,
                "new_metabolite_id": new_id,
                "name": old_name,
                "formula": old_formula,
                "old_compartment": old_comp,
                "new_compartment": new_comp,
                "affected_reaction_count": len(affected_rxns),
                "affected_reactions": " | ".join(affected_rxns),
            }
        )

    df = pd.DataFrame(records).sort_values(
        by=["action", "old_metabolite_id"],
        kind="stable",
    )
    return df, {"renamed_u": renamed, "merged_u_into_existing": merged}


def unify_quinolinate_ids_no_table(model: cobra.Model) -> Dict[str, int]:
    target_id = "quln_c"
    alias_ids = ["qns_c", "quinolinate_c", "cpd00371"]

    merged_alias = 0
    rewritten_nndpr = 0

    if target_id not in model.metabolites:
        return {"quln_alias_merged": 0, "nndpr_rewritten": 0}

    target_met = model.metabolites.get_by_id(target_id)
    for alias in alias_ids:
        if alias not in model.metabolites:
            continue
        alias_met = model.metabolites.get_by_id(alias)
        for rxn in list(alias_met.reactions):
            coeff = rxn.metabolites[alias_met]
            rxn.add_metabolites({target_met: coeff, alias_met: -coeff}, combine=True)
        model.remove_metabolites([alias_met], destructive=False)
        merged_alias += 1

    for rxn in model.reactions:
        rid = rxn.id.lower()
        rname = (rxn.name or "").lower()
        if "nndpr" not in rid and "quinolinate phosphoribosyltransferase" not in rname:
            continue
        required = ["quln_c", "prpp_c", "h_c", "nicrnt_c", "ppi_c", "co2_c"]
        if not all(mid in model.metabolites for mid in required):
            continue
        correct_stoich = {
            model.metabolites.get_by_id("quln_c"): -1.0,
            model.metabolites.get_by_id("prpp_c"): -1.0,
            model.metabolites.get_by_id("h_c"): 2.0,
            model.metabolites.get_by_id("nicrnt_c"): 1.0,
            model.metabolites.get_by_id("ppi_c"): 1.0,
            model.metabolites.get_by_id("co2_c"): 1.0,
        }
        rxn.subtract_metabolites(dict(rxn.metabolites))
        rxn.add_metabolites(correct_stoich)
        rewritten_nndpr += 1
    return {"quln_alias_merged": merged_alias, "nndpr_rewritten": rewritten_nndpr}


def run_duplicate_merge_with_audit(model: cobra.Model) -> Tuple[pd.DataFrame, Dict[str, int]]:
    grouped: Dict[Tuple[str, str, str], List[cobra.Metabolite]] = defaultdict(list)
    for met in model.metabolites:
        formula_key = canonical_formula(met.formula)
        compartment_key = get_metabolite_compartment(met)
        if formula_key and compartment_key:
            key = (compartment_key, formula_key, normalize_name(met.name or ""))
            grouped[key].append(met)

    records: List[Dict[str, object]] = []
    merge_groups = 0
    removed = 0

    for key, mets in grouped.items():
        if len(mets) <= 1:
            continue
        merge_groups += 1
        compartment, formula_key, _norm_name = key
        ordered = sorted(mets, key=lambda m: (len(m.id), m.id))
        keeper = ordered[0]
        for dup in ordered[1:]:
            affected_rxns = sorted(r.id for r in dup.reactions)
            for rxn in list(dup.reactions):
                coeff = rxn.metabolites[dup]
                rxn.add_metabolites({keeper: coeff, dup: -coeff}, combine=True)
            model.remove_metabolites([dup], destructive=False)
            removed += 1
            records.append(
                {
                    "compartment": compartment,
                    "formula": formula_key,
                    "before_metabolite_id": dup.id,
                    "after_metabolite_id": keeper.id,
                    "affected_reaction_count": len(affected_rxns),
                    "affected_reactions": " | ".join(affected_rxns),
                }
            )

    df = build_df(
        records,
        [
            "compartment",
            "formula",
            "before_metabolite_id",
            "after_metabolite_id",
            "affected_reaction_count",
            "affected_reactions",
        ],
        ["compartment", "formula", "before_metabolite_id"],
    )
    return df, {"merge_groups": merge_groups, "duplicate_nodes_removed": removed}


def add_bridge_reaction(
    model: cobra.Model,
    base_token: str,
    bridge_type: str,
    met_from: cobra.Metabolite,
    met_to: cobra.Metabolite,
    existing_ids: set[str],
) -> cobra.Reaction:
    safe_base = re.sub(r"[^A-Za-z0-9_]+", "_", base_token)
    prefix = f"TR_{safe_base}_{bridge_type.lower()}"
    rid = prefix
    idx = 1
    while rid in existing_ids or rid in model.reactions:
        idx += 1
        rid = f"{prefix}_{idx}"

    direction = {
        "cp": "c<->p",
        "ep": "e<->p",
        "ce": "c<->e",
    }.get(bridge_type.lower(), bridge_type.lower())

    rxn = cobra.Reaction(rid)
    rxn.name = f"{base_token} transport {direction}"
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000.0
    rxn.gene_reaction_rule = ""
    rxn.add_metabolites({met_from: -1.0, met_to: 1.0})
    model.add_reactions([rxn])
    existing_ids.add(rid)
    return rxn


def run_compartment_bridge_audit_with_table(model: cobra.Model) -> Tuple[pd.DataFrame, Dict[str, int]]:
    mapping: Dict[str, Dict[str, str]] = {}
    for met in model.metabolites:
        match = re.search(r"^(.*)_([cpeu])$", met.id)
        if not match:
            continue
        base, comp = match.groups()
        mapping.setdefault(base, {})[comp] = met.id

    records: List[Dict[str, object]] = []
    existing_ids = {rxn.id for rxn in model.reactions}
    cp_count = 0
    ep_count = 0
    ce_count = 0

    for base_token in sorted(mapping.keys()):
        comps = mapping[base_token]
        if "c" in comps and "p" in comps:
            m_c = model.metabolites.get_by_id(comps["c"])
            m_p = model.metabolites.get_by_id(comps["p"])
            if not set(m_c.reactions).intersection(set(m_p.reactions)):
                rxn = add_bridge_reaction(model, base_token, "cp", m_p, m_c, existing_ids)
                cp_count += 1
                records.append(
                    {
                        "reaction_id": rxn.id,
                        "base_token": base_token,
                        "bridge_type": "cp",
                        "met_from": m_p.id,
                        "met_to": m_c.id,
                        "reaction_name": rxn.name,
                        "lower_bound": rxn.lower_bound,
                        "upper_bound": rxn.upper_bound,
                        "gene_reaction_rule": rxn.gene_reaction_rule,
                        "equation": rxn.reaction,
                        "action": "",
                    }
                )

        if "p" in comps and "e" in comps:
            m_p = model.metabolites.get_by_id(comps["p"])
            m_e = model.metabolites.get_by_id(comps["e"])
            if not set(m_p.reactions).intersection(set(m_e.reactions)):
                rxn = add_bridge_reaction(model, base_token, "ep", m_e, m_p, existing_ids)
                ep_count += 1
                records.append(
                    {
                        "reaction_id": rxn.id,
                        "base_token": base_token,
                        "bridge_type": "ep",
                        "met_from": m_e.id,
                        "met_to": m_p.id,
                        "reaction_name": rxn.name,
                        "lower_bound": rxn.lower_bound,
                        "upper_bound": rxn.upper_bound,
                        "gene_reaction_rule": rxn.gene_reaction_rule,
                        "equation": rxn.reaction,
                        "action": "",
                    }
                )

        if "c" in comps and "e" in comps:
            m_c = model.metabolites.get_by_id(comps["c"])
            m_e = model.metabolites.get_by_id(comps["e"])
            if not set(m_c.reactions).intersection(set(m_e.reactions)):
                rxn = add_bridge_reaction(model, base_token, "ce", m_e, m_c, existing_ids)
                ce_count += 1
                records.append(
                    {
                        "reaction_id": rxn.id,
                        "base_token": base_token,
                        "bridge_type": "ce",
                        "met_from": m_e.id,
                        "met_to": m_c.id,
                        "reaction_name": rxn.name,
                        "lower_bound": rxn.lower_bound,
                        "upper_bound": rxn.upper_bound,
                        "gene_reaction_rule": rxn.gene_reaction_rule,
                        "equation": rxn.reaction,
                        "action": "",
                    }
                )

    df = pd.DataFrame(records).sort_values(by=["bridge_type", "base_token"], kind="stable")
    return df, {"bridge_cp": cp_count, "bridge_ep": ep_count, "bridge_ce": ce_count}


def build_df(rows: List[Dict[str, object]], columns: List[str], sort_by: List[str]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        return df
    return df.sort_values(by=sort_by, kind="stable")


def run_fba_checkpoint_table(model: cobra.Model) -> pd.DataFrame:
    solution = model.optimize()
    objective_value = 0.0
    if solution.objective_value is not None:
        objective_value = float(solution.objective_value)

    return pd.DataFrame(
        [
            {
                "status": solution.status,
                "objective_value": objective_value,
                "objective_direction": model.objective.direction,
                "objective_expression": str(model.objective.expression),
                "reaction_count": len(model.reactions),
                "metabolite_count": len(model.metabolites),
                "gene_count": len(model.genes),
            }
        ]
    )


def run_duplicate_candidate_audit_by_formula_compartment(
    model: cobra.Model,
    existing_plan_map: Dict[Tuple[str, str], str],
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    grouped: Dict[Tuple[str, str], List[cobra.Metabolite]] = defaultdict(list)
    for met in model.metabolites:
        formula_key = canonical_formula(met.formula)
        compartment_key = get_metabolite_compartment(met)
        if not formula_key or not compartment_key:
            continue
        grouped[(compartment_key, formula_key)].append(met)

    rows: List[Dict[str, object]] = []
    for (compartment_key, formula_key), mets in grouped.items():
        if len(mets) <= 1:
            continue
        sorted_mets = sorted(mets, key=lambda m: m.id)
        sorted_ids = [m.id for m in sorted_mets]
        id_name_pairs = [f"{m.id}:{(m.name or '').strip()}" for m in sorted_mets]

        default_target = sorted(sorted_mets, key=lambda m: (len(m.id), m.id))[0].id
        default_plan = f"{' | '.join(sorted_ids)} = {default_target}"
        plan_key = (compartment_key, formula_key)
        new_plan = (existing_plan_map.get(plan_key, "") or "").strip()
        # Legacy auto-filled plans are cleared so blank means "no rewrite".
        if new_plan == default_plan:
            new_plan = ""

        rows.append(
            {
                "compartment": compartment_key,
                "formula": formula_key,
                "member_count": len(sorted_ids),
                "metabolite_ids": " | ".join(sorted_ids),
                "metabolite_id_name_pairs": " | ".join(id_name_pairs),
                "new": new_plan,
            }
        )

    df = build_df(
        rows,
        [
            "compartment",
            "formula",
            "member_count",
            "metabolite_ids",
            "metabolite_id_name_pairs",
            "new",
        ],
        ["compartment", "formula"],
    )
    return df, {"duplicate_candidate_groups_formula_compartment": int(len(df))}


def load_existing_duplicate_plan_map(table_path: Path) -> Dict[Tuple[str, str], str]:
    if not table_path.exists():
        return {}
    try:
        df = pd.read_csv(table_path)
    except Exception:
        return {}

    required = {"compartment", "formula", "new"}
    if not required.issubset(df.columns):
        return {}

    plan_map: Dict[Tuple[str, str], str] = {}
    for _, row in df.iterrows():
        comp = str(row.get("compartment", "")).strip()
        formula = canonical_formula(row.get("formula", ""))
        plan = str(row.get("new", "")).strip()
        if not comp or not formula or not plan or plan.lower() == "nan":
            continue
        plan_map[(comp, formula)] = plan
    return plan_map


def parse_metabolite_id_list(text_value: object) -> List[str]:
    text = str(text_value or "")
    parts = [x.strip() for x in text.split("|")]
    out: List[str] = []
    seen: set[str] = set()
    for part in parts:
        if not part or part.lower() == "nan":
            continue
        if part not in seen:
            out.append(part)
            seen.add(part)
    return out


def parse_merge_instruction(instruction: str, default_ids: List[str]) -> Tuple[List[str], str] | None:
    raw = (instruction or "").strip()
    if not raw or raw.lower() in {"nan", "none"}:
        return None

    if "=" in raw:
        left, right = raw.split("=", 1)
        source_ids = parse_metabolite_id_list(left)
        target_id = right.strip()
        if not source_ids:
            source_ids = list(default_ids)
    else:
        target_id = raw
        source_ids = list(default_ids)

    if not target_id:
        return None

    if target_id not in source_ids:
        source_ids.append(target_id)

    ordered: List[str] = []
    seen: set[str] = set()
    for sid in source_ids:
        if sid and sid not in seen:
            ordered.append(sid)
            seen.add(sid)
    return ordered, target_id


def apply_duplicate_merge_plan_from_table(
    model: cobra.Model,
    candidates_df: pd.DataFrame,
) -> Dict[str, int]:
    planned_rows = 0
    merged_rows = 0
    merged_metabolites = 0
    skipped_rows = 0

    for _, row in candidates_df.iterrows():
        row_ids = parse_metabolite_id_list(row.get("metabolite_ids", ""))
        parsed = parse_merge_instruction(str(row.get("new", "")), row_ids)
        if parsed is None:
            continue

        planned_rows += 1
        source_ids, target_id = parsed

        if target_id not in model.metabolites:
            skipped_rows += 1
            continue

        target_met = model.metabolites.get_by_id(target_id)
        merge_sources = [sid for sid in source_ids if sid != target_id and sid in model.metabolites]
        if not merge_sources:
            continue

        merged_rows += 1
        for source_id in merge_sources:
            if source_id not in model.metabolites:
                continue
            source_met = model.metabolites.get_by_id(source_id)
            for rxn in list(source_met.reactions):
                coeff = rxn.metabolites[source_met]
                rxn.add_metabolites({target_met: coeff, source_met: -coeff}, combine=True)
            model.remove_metabolites([source_met], destructive=False)
            merged_metabolites += 1

    return {
        "merge_plan_rows": planned_rows,
        "merge_plan_rows_applied": merged_rows,
        "merged_metabolites_from_plan": merged_metabolites,
        "merge_plan_rows_skipped": skipped_rows,
    }


def load_existing_orphan_action_map(orphan_table_path: Path) -> Dict[str, str]:
    if not orphan_table_path.exists():
        return {}
    try:
        df = pd.read_csv(orphan_table_path)
    except Exception:
        return {}

    required = {"metabolite_id", "action"}
    if not required.issubset(df.columns):
        return {}

    action_map: Dict[str, str] = {}
    for _, row in df.iterrows():
        met_id = str(row.get("metabolite_id", "")).strip()
        action = str(row.get("action", "")).strip()
        if not met_id or not action or action.lower() == "nan":
            continue
        action_map[met_id] = action
    return action_map


def apply_orphan_actions_from_table(
    model: cobra.Model,
    orphan_df: pd.DataFrame,
) -> Dict[str, int]:
    orphan_action_rows = 0
    orphan_rows_deleted = 0
    orphan_action_rows_skipped = 0

    for _, row in orphan_df.iterrows():
        action = str(row.get("action", "")).strip().lower()
        if action != "no":
            continue
        orphan_action_rows += 1

        met_id = str(row.get("metabolite_id", "")).strip()
        if not met_id or met_id not in model.metabolites:
            orphan_action_rows_skipped += 1
            continue

        met = model.metabolites.get_by_id(met_id)
        model.remove_metabolites([met], destructive=False)
        orphan_rows_deleted += 1

    empty_reactions = [rxn for rxn in model.reactions if len(rxn.metabolites) == 0]
    if empty_reactions:
        model.remove_reactions(empty_reactions, remove_orphans=False)

    return {
        "orphan_action_rows": orphan_action_rows,
        "orphan_rows_deleted": orphan_rows_deleted,
        "orphan_action_rows_skipped": orphan_action_rows_skipped,
        "empty_reactions_removed_by_orphan_action": int(len(empty_reactions)),
    }


def run_post_step2_consistency_audit(
    model: cobra.Model,
    orphan_action_map: Dict[str, str] | None = None,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, int]]:
    orphan_action_map = orphan_action_map or {}

    suffix_rows: List[Dict[str, object]] = []
    mismatch_rows: List[Dict[str, object]] = []
    orphan_rows: List[Dict[str, object]] = []
    missing_formula_rows: List[Dict[str, object]] = []

    for met in model.metabolites:
        met_id = met.id
        met_formula = canonical_formula(met.formula)
        met_comp = get_metabolite_compartment(met)
        all_rxn_ids = sorted(r.id for r in met.reactions)

        if met_id.endswith("_u"):
            suffix_rows.append(
                {
                    "metabolite_id": met_id,
                    "compartment": met_comp,
                    "formula": met_formula,
                    "reaction_count": len(all_rxn_ids),
                    "reaction_ids": " | ".join(all_rxn_ids),
                }
            )

        match = re.search(r"_([A-Za-z0-9]+)$", met_id)
        id_suffix = match.group(1) if match else ""
        if id_suffix and met_comp and id_suffix != met_comp:
            mismatch_rows.append(
                {
                    "metabolite_id": met_id,
                    "id_suffix": id_suffix,
                    "compartment": met_comp,
                    "formula": met_formula,
                    "reaction_count": len(all_rxn_ids),
                }
            )

        if not met_formula:
            missing_formula_rows.append(
                {
                    "metabolite_id": met_id,
                    "compartment": met_comp,
                    "reaction_count": len(all_rxn_ids),
                }
            )

        boundary_rxns = sorted(r.id for r in met.reactions if r.boundary)
        internal_rxns = sorted(r.id for r in met.reactions if not r.boundary)

        if len(met.reactions) == 0:
            orphan_type = "no_reaction"
        elif len(internal_rxns) == 0:
            orphan_type = "boundary_only"
        elif len(internal_rxns) == 1:
            orphan_type = "single_internal_reaction"
        else:
            orphan_type = ""

        if orphan_type:
            orphan_rows.append(
                {
                    "orphan_type": orphan_type,
                    "metabolite_id": met_id,
                    "compartment": met_comp,
                    "formula": met_formula,
                    "reaction_count": len(all_rxn_ids),
                    "reaction_ids": " | ".join(all_rxn_ids),
                    "internal_reactions": " | ".join(internal_rxns),
                    "boundary_reactions": " | ".join(boundary_rxns),
                    "action": orphan_action_map.get(met_id, ""),
                }
            )

    tables = {
        "suffix_u_remaining": build_df(
            suffix_rows,
            ["metabolite_id", "compartment", "formula", "reaction_count", "reaction_ids"],
            ["metabolite_id"],
        ),
        "compartment_mismatch": build_df(
            mismatch_rows,
            ["metabolite_id", "id_suffix", "compartment", "formula", "reaction_count"],
            ["metabolite_id"],
        ),
        "orphan_metabolites": build_df(
            orphan_rows,
            [
                "orphan_type",
                "metabolite_id",
                "compartment",
                "formula",
                "reaction_count",
                "reaction_ids",
                "internal_reactions",
                "boundary_reactions",
                "action",
            ],
            ["orphan_type", "metabolite_id"],
        ),
        "missing_formula": build_df(
            missing_formula_rows,
            ["metabolite_id", "compartment", "reaction_count"],
            ["metabolite_id"],
        ),
    }

    counts = {
        "suffix_u_remaining": int(len(tables["suffix_u_remaining"])),
        "compartment_mismatch": int(len(tables["compartment_mismatch"])),
        "orphan_metabolites": int(len(tables["orphan_metabolites"])),
        "missing_formula": int(len(tables["missing_formula"])),
    }
    return tables, counts


def stringify_reactions(rxns: Iterable[str]) -> str:
    sorted_ids = sorted({x for x in rxns if x})
    return " | ".join(sorted_ids)


def make_details_string(
    reaction_ids: Iterable[str],
    reaction_map: Dict[str, Dict[str, object]],
    metabolite_names: Dict[str, str],
) -> str:
    parts: List[str] = []
    for rid in sorted({x for x in reaction_ids if x}):
        rxn = reaction_map.get(rid)
        if not rxn:
            parts.append(rid)
            continue
        name = rxn.get("name", "") if isinstance(rxn, dict) else ""
        eq = reaction_equation_from_json(rxn, metabolite_names)
        parts.append(f"{rid} [{name}] :: {eq}".strip())
    return " || ".join(parts)


def choose_best_hits(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    ordered = df.sort_values(
        by=["PID", "eVal", "COV", "bitScore"],
        ascending=[False, True, False, False],
        kind="stable",
    )
    return ordered.groupby(group_col, as_index=False).first()


def build_step1_reports(pre_step2_model: cobra.Model, pre_step2_json_path: Path) -> None:
    STEP1_DIR.mkdir(parents=True, exist_ok=True)

    pre_json = json.loads(pre_step2_json_path.read_text(encoding="utf-8"))
    ref_json = json.loads(REFERENCE_MODEL.read_text(encoding="utf-8"))
    pre_rxn_map = make_reaction_index(pre_json)
    ref_rxn_map = make_reaction_index(ref_json)
    pre_met_map = make_metabolite_name_map(pre_json)
    ref_met_map = make_metabolite_name_map(ref_json)

    gene_annotations = parse_gene_annotations(TARGET_GB)
    bbh = pd.read_csv(BBH_PARSED)
    bbh_filtered = bbh[(bbh["PID"] >= IDENTITY_THRESHOLD) & (bbh["eVal"] <= EVALUE_THRESHOLD)].copy()

    best_by_target = choose_best_hits(bbh_filtered, "subject")
    best_by_ref = choose_best_hits(bbh_filtered, "gene")

    model_gene_to_rxns = list_gene_reactions(pre_step2_model)
    ref_model = cobra.io.load_json_model(str(REFERENCE_MODEL))
    ref_gene_to_rxns = list_gene_reactions(ref_model)

    best_target_map = {row["subject"]: row for _, row in best_by_target.iterrows()}
    conflict_rows: List[Dict[str, object]] = []
    for model_gene in sorted(model_gene_to_rxns.keys()):
        hit = best_target_map.get(model_gene)
        if hit is None:
            continue
        ref_gene = str(hit["gene"])
        model_claimed = model_gene_to_rxns.get(model_gene, [])
        bigg_true = ref_gene_to_rxns.get(ref_gene, [])
        if set(model_claimed) == set(bigg_true):
            continue
        ann = gene_annotations.get(model_gene, {})
        conflict_rows.append(
            {
                "Model_Gene": model_gene,
                "Model_Gene_Symbol": ann.get("symbol", ""),
                "Model_Gene_Product": ann.get("product", ""),
                "Model_Gene_Annotation": (
                    f"symbol={ann.get('symbol', '')}; product={ann.get('product', '')}".strip("; ")
                ),
                "Best_Ref_Gene": ref_gene,
                "Identity(%)": float(hit["PID"]),
                "E_value": float(hit["eVal"]),
                "Coverage": float(hit["COV"]) if pd.notna(hit["COV"]) else "",
                "BBH": hit.get("BBH", ""),
                "Model_Claimed_Rxns": stringify_reactions(model_claimed),
                "BiGG_True_Rxns": stringify_reactions(bigg_true),
                "Model_Claimed_Rxns_Details": make_details_string(model_claimed, pre_rxn_map, pre_met_map),
                "BiGG_True_Rxns_Details": make_details_string(bigg_true, ref_rxn_map, ref_met_map),
                "Action": "",
            }
        )

    df_conflicts = pd.DataFrame(conflict_rows).sort_values(
        by=["Identity(%)", "Model_Gene"],
        ascending=[False, True],
        kind="stable",
    )
    conflicts_path = STEP1_DIR / "model_gpr_conflicts_report.csv"
    df_conflicts.to_csv(conflicts_path, index=False, encoding="utf-8-sig")

    best_ref_map = {row["gene"]: row for _, row in best_by_ref.iterrows()}
    model_gene_set = set(model_gene_to_rxns.keys())
    missing_rows: List[Dict[str, object]] = []
    for ref_gene, hit in sorted(best_ref_map.items(), key=lambda x: x[0]):
        true_rxns = ref_gene_to_rxns.get(ref_gene, [])
        if not true_rxns:
            continue
        target_gene = str(hit["subject"])
        if target_gene in model_gene_set and model_gene_to_rxns.get(target_gene):
            continue
        ann = gene_annotations.get(target_gene, {})
        reason = "target_gene_not_in_model" if target_gene not in model_gene_set else "target_gene_has_no_reactions"
        missing_rows.append(
            {
                "Ref_Gene": ref_gene,
                "Candidate_Target_Gene": target_gene,
                "Candidate_Target_Symbol": ann.get("symbol", ""),
                "Candidate_Target_Product": ann.get("product", ""),
                "Identity(%)": float(hit["PID"]),
                "E_value": float(hit["eVal"]),
                "Coverage": float(hit["COV"]) if pd.notna(hit["COV"]) else "",
                "BBH": hit.get("BBH", ""),
                "BiGG_True_Rxns": stringify_reactions(true_rxns),
                "BiGG_True_Rxns_Details": make_details_string(true_rxns, ref_rxn_map, ref_met_map),
                "Reason": reason,
                "Action": "",
            }
        )

    df_missing = pd.DataFrame(missing_rows).sort_values(
        by=["Identity(%)", "Ref_Gene"],
        ascending=[False, True],
        kind="stable",
    )
    missing_path = STEP1_DIR / "model_missing_genes_candidates.csv"
    df_missing.to_csv(missing_path, index=False, encoding="utf-8-sig")

    summary = {
        "model": str(pre_step2_json_path),
        "bbh_source": str(BBH_PARSED),
        "thresholds": {
            "identity_percent_gte": IDENTITY_THRESHOLD,
            "evalue_lte": EVALUE_THRESHOLD,
        },
        "counts": {
            "conflicts": int(len(df_conflicts)),
            "missing_candidates": int(len(df_missing)),
            "filtered_bbh_hits": int(len(bbh_filtered)),
        },
        "tables": {
            "conflicts": str(conflicts_path),
            "missing_candidates": str(missing_path),
        },
    }
    (STEP1_DIR / "step1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def run_step2_audit(pre_step2_model: cobra.Model, pre_step2_json_path: Path) -> None:
    STEP2_DIR.mkdir(parents=True, exist_ok=True)
    working = pre_step2_model.copy()

    _u_df, u_counts = run_u_suffix_cleanup_with_audit(working)
    quln_counts = unify_quinolinate_ids_no_table(working)

    duplicate_candidates_path = STEP2_DIR / "step2_table_duplicate_metabolite_candidates.csv"
    existing_plan_map = load_existing_duplicate_plan_map(duplicate_candidates_path)
    duplicate_candidates_df, candidate_counts = run_duplicate_candidate_audit_by_formula_compartment(
        working, existing_plan_map
    )
    duplicate_candidates_df.to_csv(duplicate_candidates_path, index=False, encoding="utf-8-sig")

    merge_counts = apply_duplicate_merge_plan_from_table(working, duplicate_candidates_df)
    bridge_df, bridge_counts = run_compartment_bridge_audit_with_table(working)

    orphan_path = STEP2_DIR / "step2_table_orphan_metabolites.csv"
    existing_orphan_action_map = load_existing_orphan_action_map(orphan_path)

    pre_action_tables, _ = run_post_step2_consistency_audit(working, existing_orphan_action_map)
    orphan_action_counts = apply_orphan_actions_from_table(working, pre_action_tables["orphan_metabolites"])

    post_tables, post_counts = run_post_step2_consistency_audit(working, existing_orphan_action_map)
    bridge_path = STEP2_DIR / "step2_table_compartment_bridge_reactions.csv"
    post_tables["orphan_metabolites"].to_csv(orphan_path, index=False, encoding="utf-8-sig")
    bridge_df.to_csv(bridge_path, index=False, encoding="utf-8-sig")

    cobra.io.save_json_model(working, str(FINAL_MODEL))

    for old_name in [
        "step2_table_quinolinate_related_reactions.csv",
        "step2_table_duplicate_merge_actions.csv",
        "step2_table_fba_checkpoint.csv",
        "step2_table_missing_formula.csv",
        "step2_table_suffix_u_remaining.csv",
        "step2_table_u_suffix_actions.csv",
        "step2_table_compartment_mismatch.csv",
        "DSM123_step2_duplicate_scan_model.json",
        "DSM123_step2_post_audit_model.json",
        "DSM123_step2_reconstructed_input.json",
    ]:
        old_path = STEP2_DIR / old_name
        if old_path.exists():
            old_path.unlink()

    summary = {
        "source_model_for_step2_audit": str(pre_step2_json_path),
        "updated_model_inplace": str(FINAL_MODEL),
        "notes": [
            "Duplicate candidates are grouped strictly by compartment + formula.",
            "Column 'new' drives metabolite merge actions (supports syntax like 'a | b = target').",
            "In orphan table, action=no triggers deletion; other rows are retained.",
            "Only core step2 tables are exported by request.",
        ],
        "counts": {
            **u_counts,
            **quln_counts,
            **candidate_counts,
            **merge_counts,
            **bridge_counts,
            **orphan_action_counts,
            **post_counts,
        },
        "tables": {
            "duplicate_metabolite_candidates": str(duplicate_candidates_path),
            "compartment_bridge_reactions": str(bridge_path),
            "orphan_metabolites": str(orphan_path),
        },
    }
    (STEP2_DIR / "step2_review_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> None:
    configure_cobra_cache()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    STEP1_DIR.mkdir(parents=True, exist_ok=True)
    STEP2_DIR.mkdir(parents=True, exist_ok=True)

    pre_step2_model = reconstruct_pre_step2_model()
    pre_step2_json = STEP2_DIR / "DSM123_pre_step2_reconstructed.json"
    cobra.io.save_json_model(pre_step2_model, str(pre_step2_json))

    build_step1_reports(pre_step2_model, pre_step2_json)
    run_step2_audit(pre_step2_model, pre_step2_json)

    print("Regenerated step1 + step2 review artifacts.")
    print(f"step1: {STEP1_DIR}")
    print(f"step2: {STEP2_DIR}")


if __name__ == "__main__":
    main()
