from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from .config import PipelineConfig
from .pipeline import run_pipeline


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run DSM123 metabolic reconstruction pipeline with notebook-equivalent logic.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=project_root(),
        help="Project workspace root (default: current project root).",
    )
    parser.add_argument(
        "--reference-id",
        default="Rpal_BisA53",
        help="Reference strain ID used by existing workflow files.",
    )
    parser.add_argument(
        "--target-id",
        default="consensus",
        help="Target strain ID; pipeline expects genomes/<target-id>.gb.",
    )
    parser.add_argument(
        "--target-gb",
        type=Path,
        default=None,
        help="Path to new target .gb file. If provided, it will be copied to genomes/<target-id>.gb.",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Do not clean previous outputs before running.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    cfg = PipelineConfig(
        workdir=args.workdir,
        reference_id=args.reference_id,
        target_id=args.target_id,
        target_gb=args.target_gb,
        overwrite=not args.no_overwrite,
    )

    report = run_pipeline(cfg)
    print("Pipeline completed.")
    print(f"Target GB used: {cfg.target_gb_path}")
    print(f"Final model: {cfg.final_model_path}")
    print(f"Run report: {cfg.workdir / 'pipeline_run_report.json'}")
    print(f"Final objective value: {report['fba']['final_model_objective']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
