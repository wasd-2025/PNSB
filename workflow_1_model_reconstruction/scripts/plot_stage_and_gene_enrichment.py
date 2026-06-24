from __future__ import annotations

from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DATA_DIR = ROOT / "manual_curation_outputs_merged" / "figure_data"
SVG_DIR = FIGURE_DATA_DIR / "svg"

STAGE_TABLE = FIGURE_DATA_DIR / "plot_table_stage_4d_counts.csv"
GENE_ENRICHMENT_TABLE = FIGURE_DATA_DIR / "plot_table_final_gene_enrichment_by_reaction_type.csv"

STAGE_FIG_SVG = SVG_DIR / "fig_stage_4d_grouped_bar.svg"
GENE_ENRICHMENT_FIG_SVG = SVG_DIR / "fig_final_gene_enrichment.svg"


def assert_columns(df: pd.DataFrame, required: List[str], table_path: Path) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {table_path}: {missing}")


def plot_stage_grouped_bar(stage_df: pd.DataFrame, out_svg: Path) -> None:
    ordered = stage_df.copy()
    ordered["stage"] = pd.Categorical(
        ordered["stage"],
        categories=["iDT1294Photo", "DSM123", "updated_consensus", "purple_bacteriav_DSM123"],
        ordered=True,
    )
    ordered = ordered.sort_values("stage", kind="stable")

    metrics = [
        ("genes", "Genes"),
        ("metabolites", "Metabolites"),
        ("reactions", "Reactions"),
        ("exogenous_genes_not_in_dsm123_genome", "Exogenous genes"),
    ]

    x = np.arange(len(ordered))
    width = 0.18
    offsets = np.linspace(-1.5 * width, 1.5 * width, len(metrics))
    colors = ["#2a9d8f", "#264653", "#e76f51", "#e9c46a"]

    fig, ax = plt.subplots(figsize=(12, 6.5))
    for (column, label), offset, color in zip(metrics, offsets, colors):
        values = ordered[column].astype(float).to_numpy()
        bars = ax.bar(x + offset, values, width=width, label=label, color=color, edgecolor="black", linewidth=0.5)
        for bar in bars:
            h = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                h + max(8.0, 0.01 * values.max()),
                f"{int(h)}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_title("Model Construction Dynamics Across Four Stages", fontsize=14, pad=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels([str(s) for s in ordered["stage"]], rotation=12, ha="right")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, ncol=2, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_svg, format="svg")
    plt.close(fig)


def plot_gene_enrichment(gene_df: pd.DataFrame, out_svg: Path) -> None:
    ordered = gene_df.sort_values("gene_reaction_pair_count", ascending=False, kind="stable").copy()

    x = np.arange(len(ordered))
    width = 0.36

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.5, 9.0), sharex=True, gridspec_kw={"height_ratios": [1, 1]})

    bars_pair = ax1.bar(
        x - width / 2,
        ordered["gene_reaction_pair_count"].astype(float),
        width=width,
        color="#457b9d",
        edgecolor="black",
        linewidth=0.5,
        label="Gene-reaction pairs",
    )
    bars_unique = ax1.bar(
        x + width / 2,
        ordered["unique_gene_count"].astype(float),
        width=width,
        color="#f4a261",
        edgecolor="black",
        linewidth=0.5,
        label="Unique genes",
    )

    for bars in (bars_pair, bars_unique):
        for bar in bars:
            h = bar.get_height()
            ax1.text(
                bar.get_x() + bar.get_width() / 2.0,
                h + max(2.0, 0.01 * float(ordered["gene_reaction_pair_count"].max())),
                f"{int(h)}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax1.set_title("Final Model Gene Enrichment by Reaction Type", fontsize=14, pad=10)
    ax1.set_ylabel("Count", fontsize=12)
    ax1.grid(axis="y", linestyle="--", alpha=0.35)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, loc="upper right")

    ax2.bar(
        x - width / 2,
        ordered["pair_percent"].astype(float),
        width=width,
        color="#1d3557",
        edgecolor="black",
        linewidth=0.5,
        label="Pair percent",
    )
    ax2.bar(
        x + width / 2,
        ordered["unique_gene_percent_of_model"].astype(float),
        width=width,
        color="#e63946",
        edgecolor="black",
        linewidth=0.5,
        label="Unique-gene percent of model",
    )

    for i, row in ordered.reset_index(drop=True).iterrows():
        ax2.text(i - width / 2, row["pair_percent"] + 0.8, f"{row['pair_percent']:.2f}%", ha="center", va="bottom", fontsize=9)
        ax2.text(
            i + width / 2,
            row["unique_gene_percent_of_model"] + 0.8,
            f"{row['unique_gene_percent_of_model']:.2f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax2.set_ylabel("Percent", fontsize=12)
    ax2.set_xlabel("Reaction type", fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(ordered["reaction_type"], rotation=15, ha="right")
    ax2.grid(axis="y", linestyle="--", alpha=0.35)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_svg, format="svg")
    plt.close(fig)


def main() -> None:
    SVG_DIR.mkdir(parents=True, exist_ok=True)

    stage_df = pd.read_csv(STAGE_TABLE)
    assert_columns(
        stage_df,
        ["stage", "genes", "metabolites", "reactions", "exogenous_genes_not_in_dsm123_genome"],
        STAGE_TABLE,
    )
    plot_stage_grouped_bar(stage_df, STAGE_FIG_SVG)

    gene_df = pd.read_csv(GENE_ENRICHMENT_TABLE)
    assert_columns(
        gene_df,
        [
            "reaction_type",
            "gene_reaction_pair_count",
            "pair_percent",
            "unique_gene_count",
            "unique_gene_percent_of_model",
        ],
        GENE_ENRICHMENT_TABLE,
    )
    plot_gene_enrichment(gene_df, GENE_ENRICHMENT_FIG_SVG)

    print("Generated SVG figures:")
    print(f"- {STAGE_FIG_SVG}")
    print(f"- {GENE_ENRICHMENT_FIG_SVG}")


if __name__ == "__main__":
    main()
