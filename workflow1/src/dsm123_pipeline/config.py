from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass(frozen=True)
class Thresholds:
    blast_evalue: float = 0.001
    bbh_coverage_threshold: float = 0.2
    ortho_pid_threshold: float = 65.0
    blastn_pid_threshold: float = 70.0
    blastn_aln_ratio_threshold: float = 0.8
    unannotated_pid_threshold: float = 80.0


NOTEBOOK_THRESHOLDS = Thresholds()


@dataclass
class PipelineConfig:
    workdir: Path
    reference_id: str = "Rpal_BisA53"
    target_id: str = "consensus"
    target_gb: Optional[Path] = None
    overwrite: bool = True
    reference_sbml_name: str = "iDT1294Photo.xml"
    reference_json_name: str = "iDT1294Photo.json"
    final_model_name: str = "updated_consensus.json"
    thresholds: Thresholds = NOTEBOOK_THRESHOLDS

    def __post_init__(self) -> None:
        self.workdir = Path(self.workdir).resolve()
        if self.target_gb is not None:
            self.target_gb = Path(self.target_gb).resolve()

    @property
    def genomes_dir(self) -> Path:
        return self.workdir / "genomes"

    @property
    def prots_dir(self) -> Path:
        return self.workdir / "prots"

    @property
    def nucl_dir(self) -> Path:
        return self.workdir / "nucl"

    @property
    def bbh_dir(self) -> Path:
        return self.workdir / "bbh"

    @property
    def models_dir(self) -> Path:
        return self.workdir / "Models"

    @property
    def reference_sbml_path(self) -> Path:
        return self.workdir / self.reference_sbml_name

    @property
    def reference_json_path(self) -> Path:
        return self.workdir / self.reference_json_name

    @property
    def target_gb_path(self) -> Path:
        return self.genomes_dir / f"{self.target_id}.gb"

    @property
    def target_fna_path(self) -> Path:
        return self.genomes_dir / f"{self.target_id}.fna"

    @property
    def final_model_path(self) -> Path:
        return self.workdir / self.final_model_name

    @property
    def draft_model_path(self) -> Path:
        return self.models_dir / f"{self.target_id}.json"

    def regenerated_files(self) -> List[Path]:
        return [
            self.workdir / "deletion_single_reactions(photo).csv",
            self.workdir / "flux_consensus.csv",
            self.workdir / "flux_Photo.csv",
            self.workdir / "missing_reactions.csv",
            self.workdir / "reaction_gene_relationships.csv",
            self.workdir / "updated_consensus.json",
            self.workdir / "fba_final_default_fluxes.csv",
            self.workdir / "fba_final_m9_fluxes.csv",
            self.workdir / "fba_summary.csv",
            self.workdir / "ortho_matrix.csv",
            self.workdir / "geneIDs_matrix.csv",
            self.workdir / "na_matrix.csv",
            self.workdir / "pipeline_run_report.json",
            self.draft_model_path,
            self.target_fna_path,
            self.prots_dir / f"{self.target_id}.fa",
            self.nucl_dir / f"{self.target_id}.fa",
            self.bbh_dir / f"{self.reference_id}_vs_{self.target_id}.txt",
            self.bbh_dir / f"{self.target_id}_vs_{self.reference_id}.txt",
            self.bbh_dir / f"{self.reference_id}_vs_{self.target_id}_parsed.csv",
            self.nucl_dir / f"{self.reference_id}_vs_{self.target_id}.txt",
        ]
