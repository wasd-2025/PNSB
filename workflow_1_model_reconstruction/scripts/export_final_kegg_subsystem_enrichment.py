from __future__ import annotations

import csv
import json
import re
import subprocess
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from Bio import SeqIO
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[1]
FINAL_MODEL = ROOT / "Models" / "purple_bacteriav_DSM123.json"
DSM123_GB = ROOT / "genomes" / "DSM123.gb"
CONSENSUS_GB = ROOT / "genomes" / "consensus.gb"

OUT_TABLE = ROOT / "manual_curation_outputs_merged" / "figure_data" / "plot_table_final_kegg_subsystem_enrichment.csv"
OUT_TABLE_HYBRID = (
    ROOT / "manual_curation_outputs_merged" / "figure_data" / "plot_table_final_kegg_subsystem_enrichment_hybrid.csv"
)
OUT_TABLE_HYBRID_PIE = (
    ROOT / "manual_curation_outputs_merged" / "figure_data" / "plot_table_final_kegg_subsystem_enrichment_hybrid_pie.csv"
)
OUT_SUMMARY = ROOT / "manual_curation_outputs_merged" / "figure_data" / "plot_table_final_kegg_subsystem_enrichment_summary.json"

TMP_DIR = ROOT / "manual_curation_outputs_merged" / "figure_data" / "kegg_tmp"
DSM123_FAA = TMP_DIR / "dsm123_from_gb.faa"
CONSENSUS_FAA = TMP_DIR / "consensus_from_gb.faa"
FWD_BLAST = TMP_DIR / "dsm123_vs_consensus.tsv"
REV_BLAST = TMP_DIR / "consensus_vs_dsm123.tsv"

PID_THRESHOLD = 60.0
QCOV_THRESHOLD = 0.70
SCOV_THRESHOLD = 0.70
EVALUE_THRESHOLD = 1e-5

LOGIC_WORDS = {"and", "or", "not"}
GENE_TOKEN_PATTERN = re.compile(r"[A-Za-z0-9_.:-]+")

KEGG_REST_BASE = "https://rest.kegg.jp"
ALLOWED_CLASS1 = {
    "Metabolism",
    "Environmental Information Processing",
    "Genetic Information Processing",
    "Cellular Processes",
}


def run_cmd(cmd: List[str]) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )


def parse_genes_from_gpr(gpr: str) -> List[str]:
    return [token for token in GENE_TOKEN_PATTERN.findall(gpr or "") if token.lower() not in LOGIC_WORDS]


def reaction_compartments(metabolite_ids: Iterable[str]) -> Set[str]:
    comps: Set[str] = set()
    for met_id in metabolite_ids:
        if "_" not in met_id:
            continue
        comps.add(met_id.rsplit("_", 1)[-1])
    return comps


def extract_model_genes_and_reactions(model_path: Path) -> Tuple[Set[str], List[Dict[str, object]]]:
    model_json = json.loads(model_path.read_text(encoding="utf-8"))
    genes: Set[str] = set()
    for gene in model_json.get("genes", []):
        if isinstance(gene, dict):
            gid = gene.get("id")
            if isinstance(gid, str) and gid:
                genes.add(gid)
    reactions = [r for r in model_json.get("reactions", []) if isinstance(r, dict)]
    return genes, reactions


def parse_genbank_proteins_with_ko(genbank_path: Path) -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    proteins: Dict[str, str] = {}
    ko_map: Dict[str, Set[str]] = defaultdict(set)
    for record in SeqIO.parse(str(genbank_path), "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            qualifiers = feature.qualifiers or {}
            gene_id = ""
            if qualifiers.get("locus_tag"):
                gene_id = str(qualifiers["locus_tag"][0]).strip()
            elif qualifiers.get("gene"):
                gene_id = str(qualifiers["gene"][0]).strip()
            if not gene_id:
                continue

            seq = ""
            if qualifiers.get("translation"):
                seq = str(qualifiers["translation"][0]).strip()
            else:
                try:
                    seq = str(feature.extract(record.seq).translate()).strip("*")
                except Exception:
                    seq = ""
            if seq:
                proteins[gene_id] = seq

            for dbx in qualifiers.get("db_xref", []):
                text = str(dbx).strip()
                if text.startswith("KEGG:K"):
                    ko_map[gene_id].add(text.split(":", 1)[1])
    return proteins, ko_map


def write_fasta(path: Path, seq_map: Dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for gene_id, seq in sorted(seq_map.items()):
            if not seq:
                continue
            handle.write(f">{gene_id}\n{seq}\n")


def make_blast_db(fasta_path: Path) -> None:
    run_cmd(["makeblastdb", "-in", str(fasta_path), "-dbtype", "prot"])


def run_blastp(query_faa: Path, db_faa: Path, out_tsv: Path) -> None:
    run_cmd(
        [
            "blastp",
            "-db",
            str(db_faa),
            "-query",
            str(query_faa),
            "-out",
            str(out_tsv),
            "-evalue",
            str(EVALUE_THRESHOLD),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-num_threads",
            "1",
        ]
    )


def parse_blast_best_hits(
    blast_tsv: Path,
    query_lengths: Dict[str, int],
    subject_lengths: Dict[str, int],
) -> Dict[str, Dict[str, object]]:
    best: Dict[str, Dict[str, object]] = {}
    if not blast_tsv.exists():
        return best
    with blast_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 12:
                continue
            qid, sid = parts[0], parts[1]
            try:
                pident = float(parts[2])
                aln_len = int(parts[3])
                evalue = float(parts[10])
                bitscore = float(parts[11])
            except ValueError:
                continue
            qlen = query_lengths.get(qid)
            slen = subject_lengths.get(sid)
            if not qlen or not slen:
                continue
            qcov = aln_len / qlen
            scov = aln_len / slen
            entry = {
                "query": qid,
                "subject": sid,
                "pident": pident,
                "aln_len": aln_len,
                "evalue": evalue,
                "bitscore": bitscore,
                "qcov": qcov,
                "scov": scov,
            }
            prev = best.get(qid)
            if prev is None:
                best[qid] = entry
                continue
            if (bitscore, -evalue, pident, qcov, scov) > (
                prev["bitscore"],
                -prev["evalue"],
                prev["pident"],
                prev["qcov"],
                prev["scov"],
            ):
                best[qid] = entry
    return best


def build_rbh_map(
    dsm123_proteins: Dict[str, str],
    consensus_proteins: Dict[str, str],
) -> Tuple[Dict[str, str], Dict[str, Dict[str, object]]]:
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    write_fasta(DSM123_FAA, dsm123_proteins)
    write_fasta(CONSENSUS_FAA, consensus_proteins)

    make_blast_db(CONSENSUS_FAA)
    make_blast_db(DSM123_FAA)
    run_blastp(DSM123_FAA, CONSENSUS_FAA, FWD_BLAST)
    run_blastp(CONSENSUS_FAA, DSM123_FAA, REV_BLAST)

    d_lengths = {k: len(v) for k, v in dsm123_proteins.items() if v}
    c_lengths = {k: len(v) for k, v in consensus_proteins.items() if v}

    fwd_best = parse_blast_best_hits(FWD_BLAST, d_lengths, c_lengths)
    rev_best = parse_blast_best_hits(REV_BLAST, c_lengths, d_lengths)

    rbh_map: Dict[str, str] = {}
    stats: Dict[str, Dict[str, object]] = {}
    for d_gene, hit in fwd_best.items():
        c_gene = str(hit["subject"])
        rev = rev_best.get(c_gene)
        if rev is None:
            continue
        if str(rev["subject"]) != d_gene:
            continue
        if hit["pident"] < PID_THRESHOLD:
            continue
        if hit["qcov"] < QCOV_THRESHOLD or hit["scov"] < SCOV_THRESHOLD:
            continue
        if hit["evalue"] > EVALUE_THRESHOLD:
            continue
        rbh_map[d_gene] = c_gene
        stats[d_gene] = hit
    return rbh_map, stats


def fetch_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=600) as response:
        return response.read().decode("utf-8")


def load_ko_to_pathway_map_bulk() -> Dict[str, Set[str]]:
    text = fetch_text(f"{KEGG_REST_BASE}/link/pathway/ko")
    mapping: Dict[str, Set[str]] = defaultdict(set)
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        ko_field = parts[0].strip()
        path_field = parts[1].strip()
        if not ko_field.startswith("ko:K") or not path_field.startswith("path:map"):
            continue
        ko = ko_field.split(":", 1)[1]
        path_id = path_field.split(":", 1)[1]
        # remove global/overview maps
        if path_id.startswith("map011") or path_id.startswith("map012"):
            continue
        mapping[ko].add(path_id)
    return mapping


def load_kegg_pathway_hierarchy_from_brite() -> Dict[str, Dict[str, str]]:
    text = fetch_text(f"{KEGG_REST_BASE}/get/br:br08901")
    path_info: Dict[str, Dict[str, str]] = {}
    current_class1 = ""
    current_class2 = ""
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line:
            continue
        if line.startswith("A"):
            current_class1 = line[1:].strip()
            current_class2 = ""
            continue
        if line.startswith("B"):
            current_class2 = line[1:].strip()
            continue
        if line.startswith("C"):
            payload = line[1:].strip()
            parts = payload.split()
            if not parts:
                continue
            map_code = parts[0]
            if not map_code.isdigit():
                continue
            map_id = f"map{map_code}"
            name = payload[len(map_code) :].strip()
            path_info[map_id] = {
                "name": name,
                "class1": current_class1,
                "class2": current_class2,
            }
    return path_info


def choose_majority(counter: Counter[str]) -> str:
    if not counter:
        return ""
    max_count = max(counter.values())
    ties = sorted([k for k, v in counter.items() if v == max_count])
    return ties[0] if ties else ""


def bh_fdr(pvalues: Dict[str, float]) -> Dict[str, float]:
    keys = sorted(pvalues.keys(), key=lambda k: pvalues[k])
    m = len(keys)
    if m == 0:
        return {}
    out: Dict[str, float] = {}
    prev = 1.0
    for i in range(m - 1, -1, -1):
        k = keys[i]
        rank = i + 1
        raw = pvalues[k] * m / rank
        adj = min(prev, raw, 1.0)
        out[k] = adj
        prev = adj
    return out


def classify_reaction_structure_fallback(reaction: Dict[str, object]) -> Tuple[str, str]:
    rid = str(reaction.get("id", "") or "")
    rid_upper = rid.upper()
    name = str(reaction.get("name", "") or "").lower()
    metabolites = reaction.get("metabolites", {})
    if not isinstance(metabolites, dict):
        metabolites = {}
    comps = reaction_compartments(metabolites.keys())
    if rid_upper.startswith("EX_") or rid_upper.startswith("DM_") or rid_upper.startswith("SK_"):
        return "Environmental Information Processing", "Membrane transport"
    if len(comps) > 1 or "transport" in name:
        return "Environmental Information Processing", "Membrane transport"
    return "", ""


def is_biomass_like(reaction: Dict[str, object]) -> bool:
    rid = str(reaction.get("id", "")).lower()
    name = str(reaction.get("name", "")).lower()
    if "biomass" in rid or "biomass" in name:
        return True
    obj = reaction.get("objective_coefficient", 0.0)
    try:
        return float(obj) != 0.0
    except (TypeError, ValueError):
        return False


def classify_reaction_type_fallback_detailed(reaction: Dict[str, object]) -> Tuple[str, str]:
    rid = str(reaction.get("id", "") or "")
    rid_upper = rid.upper()
    name = str(reaction.get("name", "") or "").lower()
    metabolites = reaction.get("metabolites", {})
    if not isinstance(metabolites, dict):
        metabolites = {}
    comps = reaction_compartments(metabolites.keys())

    if rid_upper.startswith("EX_"):
        return "Reaction-type fallback", "Reaction-type: Exchange"
    if rid_upper.startswith("DM_") or rid_upper.startswith("SK_"):
        return "Reaction-type fallback", "Reaction-type: Demand/Sink"
    if is_biomass_like(reaction):
        return "Reaction-type fallback", "Reaction-type: Biomass/Objective"
    if len(comps) > 1 or "transport" in name:
        return "Reaction-type fallback", "Reaction-type: Transport"
    return "Reaction-type fallback", "Reaction-type: Internal metabolic"


def main() -> None:
    model_genes, reactions = extract_model_genes_and_reactions(FINAL_MODEL)

    dsm123_proteins, _ = parse_genbank_proteins_with_ko(DSM123_GB)
    consensus_proteins, consensus_ko_by_gene = parse_genbank_proteins_with_ko(CONSENSUS_GB)

    rbh_map, _rbh_stats = build_rbh_map(dsm123_proteins, consensus_proteins)

    dsm_gene_to_kos: Dict[str, Set[str]] = defaultdict(set)
    for d_gene, c_gene in rbh_map.items():
        for ko in consensus_ko_by_gene.get(c_gene, set()):
            dsm_gene_to_kos[d_gene].add(ko)

    ko_to_paths = load_ko_to_pathway_map_bulk()
    path_info = load_kegg_pathway_hierarchy_from_brite()

    # Gene primary KEGG class assignment
    gene_to_pathways: Dict[str, Set[str]] = defaultdict(set)
    gene_to_class2: Dict[str, str] = {}
    gene_to_class1: Dict[str, str] = {}

    for gene_id, kos in dsm_gene_to_kos.items():
        class2_counter: Counter[str] = Counter()
        class1_counter: Counter[str] = Counter()
        pathway_set: Set[str] = set()
        for ko in kos:
            for pathway_id in ko_to_paths.get(ko, set()):
                info = path_info.get(pathway_id, {})
                c2 = str(info.get("class2", "")).strip()
                c1 = str(info.get("class1", "")).strip()
                if not c2 or c1 not in ALLOWED_CLASS1:
                    continue
                pathway_set.add(pathway_id)
                class2_counter[c2] += 1
                if c1:
                    class1_counter[c1] += 1
        if pathway_set:
            gene_to_pathways[gene_id] = pathway_set
        if class2_counter:
            gene_to_class2[gene_id] = choose_majority(class2_counter)
        if class1_counter:
            gene_to_class1[gene_id] = choose_majority(class1_counter)

    background_genes = sorted(gene_to_class2.keys())
    model_genes_mapped = sorted(g for g in model_genes if g in gene_to_class2)

    bg_cat_counter = Counter(gene_to_class2[g] for g in background_genes)
    model_cat_counter = Counter(gene_to_class2[g] for g in model_genes_mapped)

    N = len(background_genes)
    n = len(model_genes_mapped)

    pvals: Dict[str, float] = {}
    for cat, k in model_cat_counter.items():
        K = bg_cat_counter.get(cat, 0)
        if N == 0 or n == 0 or K == 0 or k == 0:
            pvals[cat] = 1.0
        else:
            pvals[cat] = float(hypergeom.sf(k - 1, N, K, n))
    fdr = bh_fdr(pvals)

    # Reaction assignment by mapped gene classes; fallback to membrane-transport label by structure.
    reaction_cat_counter = Counter()
    reaction_basis_counter = Counter()
    cat_to_path_counter: Dict[str, Counter[str]] = defaultdict(Counter)
    cat_basis_counter: Dict[str, Counter[str]] = defaultdict(Counter)
    cat_class1_fallback: Dict[str, str] = {}

    for rxn in reactions:
        genes = parse_genes_from_gpr(str(rxn.get("gene_reaction_rule", "") or ""))
        cat_counter: Counter[str] = Counter()
        rxn_paths: Set[str] = set()
        for g in genes:
            c2 = gene_to_class2.get(g, "")
            if c2:
                cat_counter[c2] += 1
            for p in gene_to_pathways.get(g, set()):
                rxn_paths.add(p)

        if cat_counter:
            cat = choose_majority(cat_counter)
            class1_this = choose_majority(class1_counter) if class1_counter else ""
            basis = "kegg_ko_gene_mapping"
        else:
            c1_fallback, c2_fallback = classify_reaction_type_fallback_detailed(rxn)
            cat = c2_fallback
            class1_this = c1_fallback
            basis = "reaction_type_fallback"

        reaction_cat_counter[cat] += 1
        reaction_basis_counter[basis] += 1
        cat_basis_counter[cat][basis] += 1
        if class1_this:
            cat_class1_fallback[cat] = class1_this
        for p in rxn_paths:
            info = path_info.get(p, {})
            if str(info.get("class2", "")).strip() == cat:
                cat_to_path_counter[cat][p] += 1

    total_reactions = len(reactions)

    # class2 -> class1 guess from gene assignments
    class2_to_class1_counter: Dict[str, Counter[str]] = defaultdict(Counter)
    for g, c2 in gene_to_class2.items():
        c1 = gene_to_class1.get(g, "")
        if c1:
            class2_to_class1_counter[c2][c1] += 1

    rows: List[Dict[str, object]] = []
    all_categories = sorted(set(reaction_cat_counter.keys()) | set(model_cat_counter.keys()), key=lambda c: (-reaction_cat_counter.get(c, 0), c))

    for cat in all_categories:
        rxn_count = int(reaction_cat_counter.get(cat, 0))
        rxn_pct = (rxn_count / total_reactions * 100.0) if total_reactions else 0.0

        mg = int(model_cat_counter.get(cat, 0))
        bg = int(bg_cat_counter.get(cat, 0))
        mg_pct = (mg / n * 100.0) if n else 0.0
        bg_pct = (bg / N * 100.0) if N else 0.0

        pval = pvals.get(cat, 1.0)
        qval = fdr.get(cat, 1.0)
        is_kegg_cat = cat in model_cat_counter

        class1_guess = ""
        if class2_to_class1_counter.get(cat):
            class1_guess = choose_majority(class2_to_class1_counter[cat])
        elif cat in cat_class1_fallback:
            class1_guess = cat_class1_fallback[cat]
        elif cat == "Membrane transport":
            class1_guess = "Environmental Information Processing"

        top_paths = []
        for pid, _count in cat_to_path_counter.get(cat, Counter()).most_common(5):
            info = path_info.get(pid, {})
            pname = str(info.get("name", "")).strip()
            top_paths.append(f"{pid}:{pname}" if pname else pid)

        rows.append(
            {
                "kegg_class_level1": class1_guess,
                "kegg_class_level2": cat,
                "reaction_count_final_model": rxn_count,
                "reaction_percent_final_model": round(rxn_pct, 3),
                "model_gene_count": mg,
                "model_gene_percent_of_mapped_model_genes": round(mg_pct, 3),
                "background_gene_count": bg,
                "background_gene_percent_of_mapped_genome_genes": round(bg_pct, 3),
                "enrichment_p_value": pval if is_kegg_cat else "",
                "enrichment_fdr_bh": qval if is_kegg_cat else "",
                "significant_fdr_0_05": ("yes" if qval <= 0.05 and mg > 0 else "no") if is_kegg_cat else "na",
                "reaction_assignment_basis": choose_majority(cat_basis_counter.get(cat, Counter())) if cat_basis_counter.get(cat) else "",
                "example_kegg_pathways": " | ".join(top_paths),
            }
        )

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TABLE.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "kegg_class_level1",
                "kegg_class_level2",
                "reaction_count_final_model",
                "reaction_percent_final_model",
                "model_gene_count",
                "model_gene_percent_of_mapped_model_genes",
                "background_gene_count",
                "background_gene_percent_of_mapped_genome_genes",
                "enrichment_p_value",
                "enrichment_fdr_bh",
                "significant_fdr_0_05",
                "reaction_assignment_basis",
                "example_kegg_pathways",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    # Hybrid table focused on reaction composition (for pie chart directly).
    hybrid_rows = [
        {
            "subsystem_hybrid": row["kegg_class_level2"],
            "reaction_count": row["reaction_count_final_model"],
            "percent": row["reaction_percent_final_model"],
            "assignment_basis": row["reaction_assignment_basis"],
            "class_level1": row["kegg_class_level1"],
        }
        for row in rows
    ]
    with OUT_TABLE_HYBRID.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "subsystem_hybrid",
                "reaction_count",
                "percent",
                "assignment_basis",
                "class_level1",
            ],
        )
        writer.writeheader()
        for row in hybrid_rows:
            writer.writerow(row)

    with OUT_TABLE_HYBRID_PIE.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["subsystem", "reaction_count", "percent"])
        writer.writeheader()
        for row in hybrid_rows:
            writer.writerow(
                {
                    "subsystem": row["subsystem_hybrid"],
                    "reaction_count": row["reaction_count"],
                    "percent": row["percent"],
                }
            )

    summary = {
        "final_model": str(FINAL_MODEL),
        "dsm123_genome": str(DSM123_GB),
        "consensus_genome_for_ko_transfer": str(CONSENSUS_GB),
        "rbh_thresholds": {
            "pid_gte": PID_THRESHOLD,
            "qcov_gte": QCOV_THRESHOLD,
            "scov_gte": SCOV_THRESHOLD,
            "evalue_lte": EVALUE_THRESHOLD,
        },
        "counts": {
            "final_model_genes_total": len(model_genes),
            "final_model_genes_with_kegg_class": n,
            "dsm123_genes_with_kegg_class_background": N,
            "rbh_pairs_dsm123_to_consensus": len(rbh_map),
            "consensus_genes_with_ko": int(sum(1 for _g, kos in consensus_ko_by_gene.items() if kos)),
            "unique_kos_in_background": len({ko for kos in dsm_gene_to_kos.values() for ko in kos}),
            "unique_pathways_in_background": len({p for gs in gene_to_pathways.values() for p in gs}),
            "final_model_reactions_total": total_reactions,
            "reaction_assignment_basis_counts": dict(reaction_basis_counter),
        },
        "outputs": {
            "kegg_subsystem_enrichment_table": str(OUT_TABLE),
            "kegg_subsystem_enrichment_hybrid": str(OUT_TABLE_HYBRID),
            "kegg_subsystem_enrichment_hybrid_pie": str(OUT_TABLE_HYBRID_PIE),
        },
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Exported KEGG subsystem enrichment table: {OUT_TABLE}")
    print(f"Summary: {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
