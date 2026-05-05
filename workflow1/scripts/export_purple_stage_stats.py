from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "manual_curation_outputs_merged" / "figure_data"
TARGET_GENOME_GB = ROOT / "genomes" / "DSM123.gb"

STAGE_MODELS: List[Tuple[str, Path]] = [
    ("iDT1294Photo", ROOT / "iDT1294Photo.json"),
    ("DSM123", ROOT / "Models" / "DSM123.json"),
    ("updated_consensus", ROOT / "updated_consensus.json"),
    ("purple_bacteriav_DSM123", ROOT / "Models" / "purple_bacteriav_DSM123.json"),
]

UPDATED_MODEL = ROOT / "updated_consensus.json"
FINAL_MODEL = ROOT / "Models" / "purple_bacteriav_DSM123.json"

LOGIC_WORDS = {"and", "or", "not"}
GENE_TOKEN_PATTERN = re.compile(r"[A-Za-z0-9_.:-]+")


def load_model_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: List[Dict[str, object]], columns: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def extract_gene_records(model_json: Dict[str, object]) -> Dict[str, str]:
    records: Dict[str, str] = {}
    for gene in model_json.get("genes", []):
        if not isinstance(gene, dict):
            continue
        gid = gene.get("id")
        if not isinstance(gid, str) or not gid:
            continue
        gname = str(gene.get("name", "") or "")
        if gid not in records:
            records[gid] = gname
    return records


def extract_metabolite_ids(model_json: Dict[str, object]) -> Set[str]:
    out: Set[str] = set()
    for met in model_json.get("metabolites", []):
        if not isinstance(met, dict):
            continue
        met_id = met.get("id")
        if isinstance(met_id, str) and met_id:
            out.add(met_id)
    return out


def extract_reaction_map(model_json: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        rid = rxn.get("id")
        if isinstance(rid, str) and rid:
            out[rid] = rxn
    return out


def gene_prefix(gene_id: str) -> str:
    if "_" in gene_id:
        return gene_id.split("_", 1)[0]
    return gene_id


def format_counter(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return " | ".join(f"{k}:{v}" for k, v in counter.most_common())


def load_genome_annotations(genbank_path: Path) -> Tuple[Set[str], Dict[str, Dict[str, str]]]:
    locus_tags: Set[str] = set()
    annotations: Dict[str, Dict[str, str]] = {}
    with genbank_path.open("r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qualifiers = feature.qualifiers or {}
                locus_values = qualifiers.get("locus_tag", [])
                gene_values = qualifiers.get("gene", [])
                product_values = qualifiers.get("product", [])

                if locus_values:
                    gid = str(locus_values[0]).strip()
                elif gene_values:
                    gid = str(gene_values[0]).strip()
                else:
                    gid = ""

                if not gid:
                    continue

                symbol = str(gene_values[0]).strip() if gene_values else ""
                product = str(product_values[0]).strip() if product_values else ""

                locus_tags.add(gid)
                if gid not in annotations:
                    annotations[gid] = {"symbol": symbol, "product": product}
    return locus_tags, annotations


def summarize_stage(
    stage_name: str,
    model_path: Path,
    genome_locus_tags: Set[str],
) -> Dict[str, object]:
    model_json = load_model_json(model_path)
    gene_records = extract_gene_records(model_json)
    gene_ids = sorted(gene_records.keys())

    endogenous_genes = [gid for gid in gene_ids if gid in genome_locus_tags]
    exogenous_genes = [gid for gid in gene_ids if gid not in genome_locus_tags]

    all_prefix_counts = Counter(gene_prefix(g) for g in gene_ids)
    exogenous_prefix_counts = Counter(gene_prefix(g) for g in exogenous_genes)
    endogenous_prefix_counts = Counter(gene_prefix(g) for g in endogenous_genes)

    total_genes = len(gene_ids)
    exogenous_ratio = (len(exogenous_genes) / total_genes) if total_genes else 0.0

    return {
        "stage": stage_name,
        "model_path": str(model_path),
        "genes": total_genes,
        "metabolites": len(model_json.get("metabolites", [])),
        "reactions": len(model_json.get("reactions", [])),
        "endogenous_genes_in_dsm123_genome": len(endogenous_genes),
        "exogenous_genes_not_in_dsm123_genome": len(exogenous_genes),
        "exogenous_gene_ratio": round(exogenous_ratio, 6),
        "gene_prefix_breakdown_all": format_counter(all_prefix_counts),
        "gene_prefix_breakdown_endogenous": format_counter(endogenous_prefix_counts),
        "gene_prefix_breakdown_exogenous": format_counter(exogenous_prefix_counts),
    }


def reaction_compartments(metabolite_ids: Iterable[str]) -> Set[str]:
    comps: Set[str] = set()
    for met_id in metabolite_ids:
        if "_" not in met_id:
            continue
        comp = met_id.rsplit("_", 1)[-1]
        if comp:
            comps.add(comp)
    return comps


def is_biomass_like(reaction: Dict[str, object]) -> bool:
    rid = str(reaction.get("id", "")).lower()
    rname = str(reaction.get("name", "")).lower()
    if "biomass" in rid or "biomass" in rname:
        return True
    obj = reaction.get("objective_coefficient", 0.0)
    try:
        return float(obj) != 0.0
    except (TypeError, ValueError):
        return False


def classify_reaction_type(reaction: Dict[str, object]) -> str:
    rid = str(reaction.get("id", ""))
    rid_upper = rid.upper()
    metabolites = reaction.get("metabolites", {})
    if not isinstance(metabolites, dict):
        metabolites = {}

    if rid_upper.startswith("EX_"):
        return "Exchange"
    if rid_upper.startswith("DM_") or rid_upper.startswith("SK_"):
        return "Demand/Sink"
    if is_biomass_like(reaction):
        return "Biomass/Objective"

    comps = reaction_compartments(metabolites.keys())
    if len(comps) > 1:
        return "Transport (multi-compartment)"
    if len(metabolites) <= 1:
        return "Boundary/Single-metabolite"
    return "Metabolic (internal)"


def reaction_subsystem(reaction: Dict[str, object]) -> str:
    raw = str(reaction.get("subsystem", "") or "").strip()
    return raw if raw else "Unknown"


def parse_genes_from_gpr(gpr: str) -> Set[str]:
    if not gpr:
        return set()
    out: Set[str] = set()
    for token in GENE_TOKEN_PATTERN.findall(gpr):
        if token.lower() in LOGIC_WORDS:
            continue
        out.add(token)
    return out


def normalize_gpr(gpr: str) -> str:
    return " ".join((gpr or "").split())


def summarize_reaction_types(model_json: Dict[str, object]) -> List[Dict[str, object]]:
    counts: Counter[str] = Counter()
    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        counts[classify_reaction_type(rxn)] += 1

    total = sum(counts.values())
    rows: List[Dict[str, object]] = []
    for reaction_type, count in counts.most_common():
        frac = (count / total) if total else 0.0
        rows.append(
            {
                "reaction_type": reaction_type,
                "count": count,
                "fraction": round(frac, 6),
                "percent": round(frac * 100.0, 3),
            }
        )
    return rows


def summarize_subsystems(model_json: Dict[str, object]) -> List[Dict[str, object]]:
    counts: Counter[str] = Counter()
    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        counts[reaction_subsystem(rxn)] += 1

    total = sum(counts.values())
    rows: List[Dict[str, object]] = []
    for subsystem, count in counts.most_common():
        frac = (count / total) if total else 0.0
        rows.append(
            {
                "subsystem": subsystem,
                "count": count,
                "fraction": round(frac, 6),
                "percent": round(frac * 100.0, 3),
            }
        )
    return rows


def summarize_gene_enrichment_by_reaction_groups(
    model_json: Dict[str, object],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    model_genes = set(extract_gene_records(model_json).keys())

    type_pair_counts: Counter[str] = Counter()
    type_gene_sets: Dict[str, Set[str]] = defaultdict(set)
    type_reaction_sets: Dict[str, Set[str]] = defaultdict(set)

    subsystem_pair_counts: Counter[str] = Counter()
    subsystem_gene_sets: Dict[str, Set[str]] = defaultdict(set)
    subsystem_reaction_sets: Dict[str, Set[str]] = defaultdict(set)

    total_gene_reaction_pairs = 0

    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        rid = str(rxn.get("id", "") or "")
        rtype = classify_reaction_type(rxn)
        subsystem = reaction_subsystem(rxn)
        genes = parse_genes_from_gpr(str(rxn.get("gene_reaction_rule", "") or ""))
        if not genes:
            continue

        valid_genes = sorted(g for g in genes if g in model_genes)
        if not valid_genes:
            continue

        total_gene_reaction_pairs += len(valid_genes)
        type_reaction_sets[rtype].add(rid)
        subsystem_reaction_sets[subsystem].add(rid)

        for gene_id in valid_genes:
            type_pair_counts[rtype] += 1
            type_gene_sets[rtype].add(gene_id)
            subsystem_pair_counts[subsystem] += 1
            subsystem_gene_sets[subsystem].add(gene_id)

    total_model_genes = len(model_genes)

    type_rows: List[Dict[str, object]] = []
    for rtype, pair_count in type_pair_counts.most_common():
        unique_genes = len(type_gene_sets[rtype])
        pair_frac = (pair_count / total_gene_reaction_pairs) if total_gene_reaction_pairs else 0.0
        gene_frac = (unique_genes / total_model_genes) if total_model_genes else 0.0
        type_rows.append(
            {
                "reaction_type": rtype,
                "gene_reaction_pair_count": pair_count,
                "pair_fraction": round(pair_frac, 6),
                "pair_percent": round(pair_frac * 100.0, 3),
                "unique_gene_count": unique_genes,
                "unique_gene_fraction_of_model": round(gene_frac, 6),
                "unique_gene_percent_of_model": round(gene_frac * 100.0, 3),
                "reaction_count_with_gene_rule": len(type_reaction_sets[rtype]),
            }
        )

    subsystem_rows: List[Dict[str, object]] = []
    for subsystem, pair_count in subsystem_pair_counts.most_common():
        unique_genes = len(subsystem_gene_sets[subsystem])
        pair_frac = (pair_count / total_gene_reaction_pairs) if total_gene_reaction_pairs else 0.0
        gene_frac = (unique_genes / total_model_genes) if total_model_genes else 0.0
        subsystem_rows.append(
            {
                "subsystem": subsystem,
                "gene_reaction_pair_count": pair_count,
                "pair_fraction": round(pair_frac, 6),
                "pair_percent": round(pair_frac * 100.0, 3),
                "unique_gene_count": unique_genes,
                "unique_gene_fraction_of_model": round(gene_frac, 6),
                "unique_gene_percent_of_model": round(gene_frac * 100.0, 3),
                "reaction_count_with_gene_rule": len(subsystem_reaction_sets[subsystem]),
            }
        )

    return type_rows, subsystem_rows


def build_gene_rule_index(
    model_json: Dict[str, object],
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    gene_to_reactions: Dict[str, Set[str]] = defaultdict(set)
    gene_to_gprs: Dict[str, Set[str]] = defaultdict(set)
    for rxn in model_json.get("reactions", []):
        if not isinstance(rxn, dict):
            continue
        rid = str(rxn.get("id", "") or "")
        gpr_raw = str(rxn.get("gene_reaction_rule", "") or "")
        gpr = normalize_gpr(gpr_raw)
        genes = parse_genes_from_gpr(gpr_raw)
        for gid in genes:
            gene_to_reactions[gid].add(rid)
            if gpr:
                gene_to_gprs[gid].add(gpr)
    return gene_to_reactions, gene_to_gprs


def compare_updated_to_final(
    updated_json: Dict[str, object],
    final_json: Dict[str, object],
    genome_locus_tags: Set[str],
    genome_annotations: Dict[str, Dict[str, str]],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    updated_gene_records = extract_gene_records(updated_json)
    final_gene_records = extract_gene_records(final_json)
    updated_genes = set(updated_gene_records.keys())
    final_genes = set(final_gene_records.keys())

    added_genes = sorted(final_genes - updated_genes)
    removed_genes = sorted(updated_genes - final_genes)

    final_gene_to_reactions, final_gene_to_gprs = build_gene_rule_index(final_json)

    added_gene_rows: List[Dict[str, object]] = []
    for gid in added_genes:
        rxn_ids = sorted(final_gene_to_reactions.get(gid, set()))
        gpr_rules = sorted(final_gene_to_gprs.get(gid, set()))
        ann = genome_annotations.get(gid, {})
        added_gene_rows.append(
            {
                "gene_id": gid,
                "gene_name": final_gene_records.get(gid, ""),
                "gene_symbol_from_genbank": ann.get("symbol", ""),
                "gene_product_from_genbank": ann.get("product", ""),
                "in_dsm123_genome": "yes" if gid in genome_locus_tags else "no",
                "reaction_count_in_final": len(rxn_ids),
                "reaction_ids_in_final": " | ".join(rxn_ids),
                "gpr_rules_in_final": " || ".join(gpr_rules),
            }
        )

    updated_rxn_map = extract_reaction_map(updated_json)
    final_rxn_map = extract_reaction_map(final_json)
    all_reaction_ids = sorted(set(updated_rxn_map.keys()) | set(final_rxn_map.keys()))

    gpr_change_rows: List[Dict[str, object]] = []
    gpr_changed_count = 0
    added_reaction_count = 0
    removed_reaction_count = 0

    for rid in all_reaction_ids:
        before = updated_rxn_map.get(rid)
        after = final_rxn_map.get(rid)

        if before is None and after is not None:
            change_type = "added_reaction"
            added_reaction_count += 1
            gpr_before = ""
            gpr_after = normalize_gpr(str(after.get("gene_reaction_rule", "") or ""))
            genes_before = set()
            genes_after = parse_genes_from_gpr(gpr_after)
            reaction_name = str(after.get("name", "") or "")
            subsystem = reaction_subsystem(after)
        elif after is None and before is not None:
            change_type = "removed_reaction"
            removed_reaction_count += 1
            gpr_before = normalize_gpr(str(before.get("gene_reaction_rule", "") or ""))
            gpr_after = ""
            genes_before = parse_genes_from_gpr(gpr_before)
            genes_after = set()
            reaction_name = str(before.get("name", "") or "")
            subsystem = reaction_subsystem(before)
        else:
            assert before is not None and after is not None
            gpr_before = normalize_gpr(str(before.get("gene_reaction_rule", "") or ""))
            gpr_after = normalize_gpr(str(after.get("gene_reaction_rule", "") or ""))
            if gpr_before == gpr_after:
                continue
            change_type = "gpr_changed"
            gpr_changed_count += 1
            genes_before = parse_genes_from_gpr(gpr_before)
            genes_after = parse_genes_from_gpr(gpr_after)
            reaction_name = str(after.get("name", "") or before.get("name", "") or "")
            subsystem = reaction_subsystem(after if after is not None else before)

        added_genes_in_rule = sorted(genes_after - genes_before)
        removed_genes_in_rule = sorted(genes_before - genes_after)

        gpr_change_rows.append(
            {
                "change_type": change_type,
                "reaction_id": rid,
                "reaction_name": reaction_name,
                "subsystem": subsystem,
                "gpr_before": gpr_before,
                "gpr_after": gpr_after,
                "genes_before": " | ".join(sorted(genes_before)),
                "genes_after": " | ".join(sorted(genes_after)),
                "added_genes_in_rule": " | ".join(added_genes_in_rule),
                "removed_genes_in_rule": " | ".join(removed_genes_in_rule),
            }
        )

    gpr_change_rows.sort(key=lambda x: (x["change_type"], x["reaction_id"]))

    updated_met_ids = extract_metabolite_ids(updated_json)
    final_met_ids = extract_metabolite_ids(final_json)

    summary_rows = [
        {
            "updated_model_path": str(UPDATED_MODEL),
            "final_model_path": str(FINAL_MODEL),
            "genes_updated": len(updated_genes),
            "genes_final": len(final_genes),
            "genes_added": len(added_genes),
            "genes_removed": len(removed_genes),
            "reactions_updated": len(updated_rxn_map),
            "reactions_final": len(final_rxn_map),
            "reactions_added": len(set(final_rxn_map) - set(updated_rxn_map)),
            "reactions_removed": len(set(updated_rxn_map) - set(final_rxn_map)),
            "gpr_changed_reactions": gpr_changed_count,
            "added_reaction_rows_in_gpr_table": added_reaction_count,
            "removed_reaction_rows_in_gpr_table": removed_reaction_count,
            "metabolites_updated": len(updated_met_ids),
            "metabolites_final": len(final_met_ids),
            "metabolites_added": len(final_met_ids - updated_met_ids),
            "metabolites_removed": len(updated_met_ids - final_met_ids),
            "added_genes_with_reaction_links": sum(
                1 for gid in added_genes if len(final_gene_to_reactions.get(gid, set())) > 0
            ),
        }
    ]
    return added_gene_rows, gpr_change_rows, summary_rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    genome_locus_tags, genome_annotations = load_genome_annotations(TARGET_GENOME_GB)

    stage_rows: List[Dict[str, object]] = []
    stage_model_paths: Dict[str, str] = {}
    for stage_name, model_path in STAGE_MODELS:
        if not model_path.exists():
            continue
        stage_rows.append(summarize_stage(stage_name, model_path, genome_locus_tags))
        stage_model_paths[stage_name] = str(model_path)

    stage_table = OUT_DIR / "plot_table_stage_4d_counts.csv"
    write_csv(
        stage_table,
        stage_rows,
        [
            "stage",
            "model_path",
            "genes",
            "metabolites",
            "reactions",
            "endogenous_genes_in_dsm123_genome",
            "exogenous_genes_not_in_dsm123_genome",
            "exogenous_gene_ratio",
            "gene_prefix_breakdown_all",
            "gene_prefix_breakdown_endogenous",
            "gene_prefix_breakdown_exogenous",
        ],
    )

    final_json = load_model_json(FINAL_MODEL)
    updated_json = load_model_json(UPDATED_MODEL)

    reaction_type_rows = summarize_reaction_types(final_json)
    reaction_type_table = OUT_DIR / "plot_table_final_reaction_type_enrichment.csv"
    write_csv(
        reaction_type_table,
        reaction_type_rows,
        ["reaction_type", "count", "fraction", "percent"],
    )

    subsystem_rows = summarize_subsystems(final_json)
    subsystem_table = OUT_DIR / "plot_table_final_reaction_subsystem_enrichment.csv"
    write_csv(
        subsystem_table,
        subsystem_rows,
        ["subsystem", "count", "fraction", "percent"],
    )

    gene_type_rows, gene_subsystem_rows = summarize_gene_enrichment_by_reaction_groups(final_json)
    gene_type_table = OUT_DIR / "plot_table_final_gene_enrichment_by_reaction_type.csv"
    gene_subsystem_table = OUT_DIR / "plot_table_final_gene_enrichment_by_subsystem.csv"
    write_csv(
        gene_type_table,
        gene_type_rows,
        [
            "reaction_type",
            "gene_reaction_pair_count",
            "pair_fraction",
            "pair_percent",
            "unique_gene_count",
            "unique_gene_fraction_of_model",
            "unique_gene_percent_of_model",
            "reaction_count_with_gene_rule",
        ],
    )
    write_csv(
        gene_subsystem_table,
        gene_subsystem_rows,
        [
            "subsystem",
            "gene_reaction_pair_count",
            "pair_fraction",
            "pair_percent",
            "unique_gene_count",
            "unique_gene_fraction_of_model",
            "unique_gene_percent_of_model",
            "reaction_count_with_gene_rule",
        ],
    )

    added_gene_rows, gpr_change_rows, delta_summary_rows = compare_updated_to_final(
        updated_json, final_json, genome_locus_tags, genome_annotations
    )
    added_genes_table = OUT_DIR / "plot_table_updated_to_final_added_genes_with_gpr.csv"
    gpr_changes_table = OUT_DIR / "plot_table_updated_to_final_gpr_changes.csv"
    delta_summary_table = OUT_DIR / "plot_table_updated_to_final_overall_delta.csv"

    write_csv(
        added_genes_table,
        added_gene_rows,
        [
            "gene_id",
            "gene_name",
            "gene_symbol_from_genbank",
            "gene_product_from_genbank",
            "in_dsm123_genome",
            "reaction_count_in_final",
            "reaction_ids_in_final",
            "gpr_rules_in_final",
        ],
    )
    write_csv(
        gpr_changes_table,
        gpr_change_rows,
        [
            "change_type",
            "reaction_id",
            "reaction_name",
            "subsystem",
            "gpr_before",
            "gpr_after",
            "genes_before",
            "genes_after",
            "added_genes_in_rule",
            "removed_genes_in_rule",
        ],
    )
    write_csv(
        delta_summary_table,
        delta_summary_rows,
        list(delta_summary_rows[0].keys()),
    )

    summary = {
        "dsm123_genome_for_endogenous_exogenous_definition": str(TARGET_GENOME_GB),
        "stage_models": stage_model_paths,
        "outputs": {
            "stage_4d_counts": str(stage_table),
            "final_reaction_type_enrichment": str(reaction_type_table),
            "final_reaction_subsystem_enrichment": str(subsystem_table),
            "final_gene_enrichment_by_reaction_type": str(gene_type_table),
            "final_gene_enrichment_by_subsystem": str(gene_subsystem_table),
            "updated_to_final_added_genes_with_gpr": str(added_genes_table),
            "updated_to_final_gpr_changes": str(gpr_changes_table),
            "updated_to_final_overall_delta": str(delta_summary_table),
        },
        "counts": {
            "stage_rows": len(stage_rows),
            "final_reaction_type_groups": len(reaction_type_rows),
            "final_subsystem_groups": len(subsystem_rows),
            "final_gene_reaction_type_groups": len(gene_type_rows),
            "final_gene_subsystem_groups": len(gene_subsystem_rows),
            "updated_to_final_added_gene_rows": len(added_gene_rows),
            "updated_to_final_gpr_change_rows": len(gpr_change_rows),
        },
    }
    (OUT_DIR / "plot_table_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Exported figure data tables:")
    for key, value in summary["outputs"].items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()
