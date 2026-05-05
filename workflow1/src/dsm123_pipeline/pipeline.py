from __future__ import annotations

import csv
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import appdirs
import pandas as pd
from Bio import SeqIO

from .config import PipelineConfig


BLAST_COLS = [
    "gene",
    "subject",
    "PID",
    "alnLength",
    "mismatchCount",
    "gapOpenCount",
    "queryStart",
    "queryEnd",
    "subjectStart",
    "subjectEnd",
    "eVal",
    "bitScore",
]


def configure_cobra_cache(workdir: Path) -> None:
    cache_dir = workdir / ".cobra_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    appdirs.user_cache_dir = lambda *args, **kwargs: str(cache_dir)


def run_cmd(cmd: List[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def ensure_dirs(cfg: PipelineConfig) -> None:
    for folder in [cfg.genomes_dir, cfg.prots_dir, cfg.nucl_dir, cfg.bbh_dir, cfg.models_dir]:
        folder.mkdir(parents=True, exist_ok=True)


def clean_outputs(cfg: PipelineConfig) -> None:
    for path in cfg.regenerated_files():
        if path.exists() and path.is_file():
            path.unlink()

    for artifact in cfg.prots_dir.glob(f"{cfg.target_id}.fa.*"):
        if artifact.is_file():
            artifact.unlink()

    for artifact in cfg.genomes_dir.glob(f"{cfg.target_id}.fna.*"):
        if artifact.is_file():
            artifact.unlink()


def stage_target_genome(cfg: PipelineConfig) -> None:
    if cfg.target_gb is not None:
        if not cfg.target_gb.exists():
            raise FileNotFoundError(f"Target gb file not found: {cfg.target_gb}")
        shutil.copy2(cfg.target_gb, cfg.target_gb_path)

    if not cfg.target_gb_path.exists():
        raise FileNotFoundError(
            f"Target gb file not found at {cfg.target_gb_path}. "
            f"Please place your file there or pass --target-gb."
        )


def parse_genome(genbank_file: Path, output_fasta: Path, nucleotide: bool) -> None:
    count = 0
    with genbank_file.open("r", encoding="utf-8") as handle, output_fasta.open("w", encoding="utf-8") as fout:
        for record in SeqIO.parse(handle, "genbank"):
            fallback_index = 0
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                try:
                    seq = feature.extract(record.seq)
                    seq_text = str(seq) if nucleotide else str(seq.translate())
                except Exception:
                    continue

                if "locus_tag" in feature.qualifiers:
                    gene_id = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers:
                    gene_id = feature.qualifiers["gene"][0]
                else:
                    gene_id = f"gene_{fallback_index}"
                    fallback_index += 1
                fout.write(f">{gene_id}\n{seq_text}\n")
                count += 1

    if count == 0:
        raise RuntimeError(f"No CDS exported from {genbank_file}")


def ensure_reference_fastas(cfg: PipelineConfig) -> None:
    ref_prot = cfg.prots_dir / f"{cfg.reference_id}.fa"
    ref_nucl = cfg.nucl_dir / f"{cfg.reference_id}.fa"

    if ref_prot.exists() and ref_prot.stat().st_size > 0 and ref_nucl.exists() and ref_nucl.stat().st_size > 0:
        return

    reference_gb = cfg.genomes_dir / f"{cfg.reference_id}.gb"
    if reference_gb.exists():
        try:
            parse_genome(reference_gb, ref_nucl, nucleotide=True)
        except RuntimeError:
            pass
        try:
            parse_genome(reference_gb, ref_prot, nucleotide=False)
        except RuntimeError:
            pass

    if ref_prot.exists() and ref_prot.stat().st_size > 0 and ref_nucl.exists() and ref_nucl.stat().st_size > 0:
        return

    fallback_dirs = [
        cfg.workdir.parent / "workflow" / "data" / "raw",
        cfg.workdir.parent.parent / "Rhodopseudomonas palustris models" / "dsm123" / "workflow",
        cfg.workdir.parent.parent / "Rhodopseudomonas palustris models" / "dsm130" / "workflow",
    ]
    for root in fallback_dirs:
        src_prot = root / "prots" / f"{cfg.reference_id}.fa"
        src_nucl = root / "nucl" / f"{cfg.reference_id}.fa"
        if src_prot.exists() and src_prot.stat().st_size > 0 and not (ref_prot.exists() and ref_prot.stat().st_size > 0):
            shutil.copy2(src_prot, ref_prot)
        if src_nucl.exists() and src_nucl.stat().st_size > 0 and not (ref_nucl.exists() and ref_nucl.stat().st_size > 0):
            shutil.copy2(src_nucl, ref_nucl)

    if not (ref_prot.exists() and ref_prot.stat().st_size > 0 and ref_nucl.exists() and ref_nucl.stat().st_size > 0):
        raise FileNotFoundError("Reference FASTA files are missing or empty.")


def make_blast_db(input_file: Path, db_type: str) -> None:
    run_cmd(["makeblastdb", "-in", str(input_file), "-dbtype", db_type])


def run_blastp(query_fa: Path, db_fa: Path, output_txt: Path, evalue: float) -> None:
    run_cmd(
        [
            "blastp",
            "-db",
            str(db_fa),
            "-query",
            str(query_fa),
            "-out",
            str(output_txt),
            "-evalue",
            str(evalue),
            "-outfmt",
            "6",
            "-num_threads",
            "1",
        ]
    )


def run_blastn(query_fa: Path, db_fna: Path, output_txt: Path, evalue: float) -> None:
    run_cmd(
        [
            "blastn",
            "-db",
            str(db_fna),
            "-query",
            str(query_fa),
            "-out",
            str(output_txt),
            "-evalue",
            str(evalue),
            "-outfmt",
            "6",
            "-num_threads",
            "1",
        ]
    )


def get_gene_lengths(fasta_file: Path) -> pd.DataFrame:
    with fasta_file.open("r", encoding="utf-8") as handle:
        return pd.DataFrame([{"gene": r.name, "gene_length": len(r.seq)} for r in SeqIO.parse(handle, "fasta")])


def build_bbh_table(
    forward_txt: Path,
    reverse_txt: Path,
    query_fa: Path,
    subject_fa: Path,
    coverage_threshold: float,
) -> pd.DataFrame:
    query_lengths = get_gene_lengths(query_fa)
    subject_lengths = get_gene_lengths(subject_fa)

    bbh = pd.read_csv(forward_txt, sep="\t", names=BLAST_COLS)
    bbh = pd.merge(bbh, query_lengths, on="gene", how="left")
    bbh["COV"] = bbh["alnLength"] / bbh["gene_length"]
    bbh = bbh[bbh["COV"] >= coverage_threshold]

    bbh2 = pd.read_csv(reverse_txt, sep="\t", names=BLAST_COLS)
    bbh2 = pd.merge(bbh2, subject_lengths, on="gene", how="left")
    bbh2["COV"] = bbh2["alnLength"] / bbh2["gene_length"]
    bbh2 = bbh2[bbh2["COV"] >= coverage_threshold]

    out = pd.DataFrame()
    for gene in bbh["gene"].unique():
        res = bbh[bbh["gene"] == gene]
        if res.empty:
            continue
        best_hit = res.loc[res["PID"].idxmax()].copy()
        best_gene = best_hit["subject"]

        res2 = bbh2[bbh2["gene"] == best_gene]
        if res2.empty:
            continue
        best_hit2 = res2.loc[res2["PID"].idxmax()]
        best_hit["BBH"] = "<=>" if gene == best_hit2["subject"] else "->"
        out = pd.concat([out, pd.DataFrame(best_hit).transpose()], ignore_index=True)
    return out


def build_orthology_matrices(parsed_bbh: pd.DataFrame, cfg: PipelineConfig) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    from cobra.io import load_json_model

    model = load_json_model(str(cfg.reference_json_path))
    reference_genes = [g.id for g in model.genes]

    ortho_pid = pd.DataFrame(index=reference_genes, columns=[cfg.target_id])
    gene_ids = pd.DataFrame(index=reference_genes, columns=[cfg.target_id])

    pid_values: List[float] = []
    id_values: List[str] = []
    for gene in reference_genes:
        current = parsed_bbh[parsed_bbh["gene"] == gene].reset_index(drop=True)
        if current.empty:
            pid_values.append(0.0)
            id_values.append("None")
        else:
            pid_values.append(float(current.iloc[0]["PID"]))
            id_values.append(str(current.iloc[0]["subject"]))

    ortho_pid[cfg.target_id] = pid_values
    gene_ids[cfg.target_id] = id_values

    ortho_bin = (ortho_pid > cfg.thresholds.ortho_pid_threshold).astype(int)
    return ortho_bin, gene_ids, reference_genes


def gbk_to_fna(gbk_file: Path, out_fna: Path) -> None:
    with gbk_file.open("r", encoding="utf-8") as input_handle, out_fna.open("w", encoding="utf-8") as output_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            output_handle.write(f">{record.id} {record.description}\n{record.seq}\n")


def extract_seq(fasta_file: Path, contig: str, start: int, end: int) -> str:
    with fasta_file.open("r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.name != contig:
                continue
            if end > start:
                section = record[start:end]
            else:
                section = record[end - 1 : start + 1].reverse_complement()
            return str(section.seq)
    return ""


def add_unannotated_orfs(
    ortho_matrix: pd.DataFrame,
    gene_ids_matrix: pd.DataFrame,
    na_matrix: pd.DataFrame,
    reference_genes: List[str],
    cfg: PipelineConfig,
) -> Dict[str, str]:
    if "contig_1" not in na_matrix.columns:
        return {}

    na_model_genes = na_matrix.drop([g for g in na_matrix.index if g not in reference_genes])
    present_genes = ortho_matrix[cfg.target_id][ortho_matrix[cfg.target_id] == 1].index.tolist()
    candidate_genes = na_model_genes["contig_1"][
        na_model_genes["contig_1"] >= cfg.thresholds.unannotated_pid_threshold
    ].index.tolist()
    unannotated = set(candidate_genes) - set(present_genes)

    blastn_file = cfg.nucl_dir / f"{cfg.reference_id}_vs_{cfg.target_id}.txt"
    blastn_df = pd.read_csv(blastn_file, sep="\t", names=BLAST_COLS)
    pseudogenes: Dict[str, str] = {}

    for _, row in blastn_df[blastn_df["gene"].isin(unannotated)].iterrows():
        gene = row["gene"]
        seq = extract_seq(
            cfg.target_fna_path,
            str(row["subject"]),
            int(row["subjectStart"]) - 1,
            int(row["subjectEnd"]),
        )
        if "*" in seq:
            pseudogenes[gene] = seq
            if gene not in ortho_matrix.index:
                ortho_matrix.loc[gene, cfg.target_id] = 1
                gene_ids_matrix.loc[gene, cfg.target_id] = f"{gene}_ortholog"

    return pseudogenes


def build_draft_model(ortho_matrix: pd.DataFrame, gene_ids_matrix: pd.DataFrame, cfg: PipelineConfig) -> Path:
    import cobra
    from cobra.io import load_json_model
    from cobra.manipulation.delete import remove_genes
    from cobra.manipulation.modify import rename_genes

    ref_model = load_json_model(str(cfg.reference_json_path))
    non_homologous = ortho_matrix[cfg.target_id][ortho_matrix[cfg.target_id] == 0].index.tolist()

    to_delete = []
    for gene in non_homologous:
        try:
            to_delete.append(ref_model.genes.get_by_id(gene))
        except KeyError:
            continue

    draft = ref_model.copy()
    remove_genes(draft, to_delete, remove_reactions=True)
    draft.id = cfg.target_id

    mapping = {
        str(old): str(new)
        for old, new in gene_ids_matrix[cfg.target_id].to_dict().items()
        if pd.notna(new) and str(new) != "None"
    }
    rename_genes(draft, mapping)

    cobra.io.save_json_model(draft, str(cfg.draft_model_path))
    return cfg.draft_model_path


def gapfill_with_flux_filter(cfg: PipelineConfig, draft_path: Path) -> Tuple[int, int, float]:
    import cobra
    from cobra.flux_analysis import single_reaction_deletion
    from cobra.io import load_json_model

    ref_model = load_json_model(str(cfg.reference_json_path))
    draft_model = load_json_model(str(draft_path))

    single_reaction_deletion(ref_model, processes=1).to_csv(cfg.workdir / "deletion_single_reactions(photo).csv")

    draft_model.optimize().to_frame().to_csv(cfg.workdir / "flux_consensus.csv")
    ref_model.optimize().to_frame().to_csv(cfg.workdir / "flux_Photo.csv")

    ref_flux = pd.read_csv(cfg.workdir / "flux_Photo.csv")
    ref_flux.columns = ["reactions" if i == 0 else c for i, c in enumerate(ref_flux.columns)]
    filtered = ref_flux[(ref_flux["fluxes"] != 0) | (ref_flux["reduced_costs"] != 0)]

    missing = set(filtered["reactions"]) - {r.id for r in draft_model.reactions}
    filtered[filtered["reactions"].isin(missing)][["reactions", "fluxes", "reduced_costs"]].to_csv(
        cfg.workdir / "missing_reactions.csv",
        index=False,
    )

    added: Dict[str, List[str]] = {}
    for rid in sorted(missing):
        try:
            rxn = ref_model.reactions.get_by_id(rid).copy()
        except KeyError:
            continue
        draft_model.add_reactions([rxn])
        added[rid] = [g.id for g in rxn.genes]

    cobra.io.save_json_model(draft_model, str(cfg.final_model_path))

    with (cfg.workdir / "reaction_gene_relationships.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Reaction", "Genes"])
        writer.writeheader()
        for rid, genes in added.items():
            writer.writerow({"Reaction": rid, "Genes": ",".join(genes)})

    final_objective = load_json_model(str(cfg.final_model_path)).optimize().objective_value
    return len(missing), len(added), float(final_objective)


def apply_m9(model) -> None:
    for reaction in model.reactions:
        if reaction.id.startswith("EX_"):
            reaction.lower_bound = 0

    fixed_lb = {
        "EX_mobd_e": -1000,
        "EX_cu2_e": -1000,
        "EX_bo3_e": -1000,
        "EX_na1_e": -1000,
        "EX_mn2_e": -1000,
        "EX_zn2_e": -1000,
        "EX_k_e": -1000,
        "EX_cobalt2_e": -1000,
        "EX_ca2_e": -1000,
        "EX_cl_e": -1000,
        "EX_mg2_e": -1000,
        "EX_fe3_e": -10,
        "EX_so4_e": -1000,
        "EX_for_e": -10,
        "EX_nh4_e": -5,
        "EX_o2_e": -20,
        "EX_n2_e": -2,
    }
    for rid, lb in fixed_lb.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = lb


def run_simple_fba(cfg: PipelineConfig) -> Tuple[float, float]:
    from cobra.io import load_json_model

    final_model = load_json_model(str(cfg.final_model_path))
    sol_default = final_model.optimize()
    sol_default.to_frame().to_csv(cfg.workdir / "fba_final_default_fluxes.csv")

    final_m9 = load_json_model(str(cfg.final_model_path))
    apply_m9(final_m9)
    sol_m9 = final_m9.optimize()
    sol_m9.to_frame().to_csv(cfg.workdir / "fba_final_m9_fluxes.csv")

    pd.DataFrame(
        [
            {"condition": "default", "objective_value": sol_default.objective_value},
            {"condition": "m9", "objective_value": sol_m9.objective_value},
        ]
    ).to_csv(cfg.workdir / "fba_summary.csv", index=False)

    return float(sol_default.objective_value), float(sol_m9.objective_value)


def run_pipeline(cfg: PipelineConfig) -> Dict[str, object]:
    configure_cobra_cache(cfg.workdir)
    from cobra.io import read_sbml_model, save_json_model
    ensure_dirs(cfg)

    if cfg.overwrite:
        clean_outputs(cfg)

    stage_target_genome(cfg)

    if not cfg.reference_json_path.exists():
        save_json_model(read_sbml_model(str(cfg.reference_sbml_path)), str(cfg.reference_json_path))

    parse_genome(cfg.target_gb_path, cfg.prots_dir / f"{cfg.target_id}.fa", nucleotide=False)
    parse_genome(cfg.target_gb_path, cfg.nucl_dir / f"{cfg.target_id}.fa", nucleotide=True)

    ensure_reference_fastas(cfg)

    ref_prot = cfg.prots_dir / f"{cfg.reference_id}.fa"
    target_prot = cfg.prots_dir / f"{cfg.target_id}.fa"
    make_blast_db(target_prot, "prot")
    make_blast_db(ref_prot, "prot")

    ref_vs_target = cfg.bbh_dir / f"{cfg.reference_id}_vs_{cfg.target_id}.txt"
    target_vs_ref = cfg.bbh_dir / f"{cfg.target_id}_vs_{cfg.reference_id}.txt"
    run_blastp(ref_prot, target_prot, ref_vs_target, cfg.thresholds.blast_evalue)
    run_blastp(target_prot, ref_prot, target_vs_ref, cfg.thresholds.blast_evalue)

    parsed_bbh = build_bbh_table(
        ref_vs_target,
        target_vs_ref,
        ref_prot,
        target_prot,
        cfg.thresholds.bbh_coverage_threshold,
    )
    parsed_bbh.to_csv(cfg.bbh_dir / f"{cfg.reference_id}_vs_{cfg.target_id}_parsed.csv", index=False)

    ortho_matrix, gene_ids_matrix, reference_genes = build_orthology_matrices(parsed_bbh, cfg)

    gbk_to_fna(cfg.target_gb_path, cfg.target_fna_path)
    make_blast_db(cfg.target_fna_path, "nucl")

    ref_nucl = cfg.nucl_dir / f"{cfg.reference_id}.fa"
    ref_vs_target_nucl = cfg.nucl_dir / f"{cfg.reference_id}_vs_{cfg.target_id}.txt"
    run_blastn(ref_nucl, cfg.target_fna_path, ref_vs_target_nucl, cfg.thresholds.blast_evalue)

    na_df = pd.read_csv(ref_vs_target_nucl, sep="\t", names=BLAST_COLS)
    na_df = na_df[
        (na_df["PID"] > cfg.thresholds.blastn_pid_threshold)
        & (na_df["alnLength"] > cfg.thresholds.blastn_aln_ratio_threshold * na_df["queryEnd"])
    ]
    na_df = na_df.groupby("gene").first().reset_index()
    na_df.to_csv(cfg.workdir / "na_matrix.csv", index=False)
    na_matrix = pd.pivot_table(na_df, index="gene", columns="subject", values="PID")

    pseudogenes = add_unannotated_orfs(ortho_matrix, gene_ids_matrix, na_matrix, reference_genes, cfg)

    ortho_matrix.to_csv(cfg.workdir / "ortho_matrix.csv")
    gene_ids_matrix.to_csv(cfg.workdir / "geneIDs_matrix.csv")

    draft_path = build_draft_model(ortho_matrix, gene_ids_matrix, cfg)
    missing_count, added_count, final_obj = gapfill_with_flux_filter(cfg, draft_path)
    default_obj, m9_obj = run_simple_fba(cfg)

    report = {
        "run_timestamp": datetime.now().isoformat(timespec="seconds"),
        "reference_id": cfg.reference_id,
        "target_id": cfg.target_id,
        "target_gb_path": str(cfg.target_gb_path),
        "thresholds": {
            "bbh_cov_threshold": cfg.thresholds.bbh_coverage_threshold,
            "ortho_pid_threshold": cfg.thresholds.ortho_pid_threshold,
            "blastn_pid_threshold": cfg.thresholds.blastn_pid_threshold,
            "blastn_alignment_ratio_threshold": cfg.thresholds.blastn_aln_ratio_threshold,
            "unannotated_pid_threshold": cfg.thresholds.unannotated_pid_threshold,
            "blast_evalue": cfg.thresholds.blast_evalue,
        },
        "gapfill": {
            "missing_reactions": missing_count,
            "added_reactions": added_count,
        },
        "fba": {
            "default_objective": default_obj,
            "m9_objective": m9_obj,
            "final_model_objective": final_obj,
        },
        "pseudogene_hits": pseudogenes,
    }

    (cfg.workdir / "pipeline_run_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


