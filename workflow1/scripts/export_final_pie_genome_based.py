from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
FINAL_MODEL = ROOT / "Models" / "purple_bacteriav_DSM123.json"
GENOME_GB = ROOT / "genomes" / "DSM123.gb"
OUT_TABLE = ROOT / "manual_curation_outputs_merged" / "figure_data" / "plot_table_final_pie_genome_based_subsystem.csv"

LOGIC_WORDS = {"and", "or", "not"}
GENE_TOKEN_PATTERN = re.compile(r"[A-Za-z0-9_.:-]+")

CATEGORY_ORDER = [
    "Pentose Phosphate Pathway",
    "Nitrogen Metabolism",
    "Carbohydrate Metabolism",
    "Pyruvate Metabolism",
    "Energy Metabolism",
    "Carotenoid Biosynthesis",
    "Aromatic Compounds Metabolism",
    "Photosynthesis",
    "Nucleotide Metabolism",
    "Alternate Carbon Metabolism",
    "Cell Envelope Biosynthesis",
    "Amino Acid Metabolism",
    "Extracellular exchange",
    "Cofactor and Prosthetic Group Biosynthesis",
    "Lipid Metabolism",
    "Transport",
    "Others",
]

KEYWORDS: Dict[str, List[str]] = {
    "Photosynthesis": [
        "photosynth",
        "reaction center",
        "light harvesting",
        "bacteriochlorophyll",
        "chlorophyll",
        "photon",
        "puf",
        "puh",
        "psb",
        "psi",
        "photochemical",
    ],
    "Carotenoid Biosynthesis": [
        "carotenoid",
        "carotene",
        "lycopene",
        "phytoene",
        "crt",
    ],
    "Aromatic Compounds Metabolism": [
        "aromatic",
        "benzoate",
        "catechol",
        "protocatechu",
        "phenylacet",
        "phenylprop",
        "shikimate",
        "chorismate",
    ],
    "Pentose Phosphate Pathway": [
        "pentose phosphate",
        "transketolase",
        "transaldolase",
        "ribulose 5",
        "ribose 5",
        "6-phosphogluconate",
        "phosphogluconate",
        "gluconolactone",
        "zwf",
        "gnd",
        "rpe",
        "rpi",
    ],
    "Pyruvate Metabolism": [
        "pyruvate",
        "pyruv",
    ],
    "Nitrogen Metabolism": [
        "nitrogen",
        "nitrate",
        "nitrite",
        "ammonia",
        "ammonium",
        "urea",
        "urease",
        "cyanate",
        "nitrile",
        "nitric",
        "dinitrogen",
        "glutamine synthetase",
    ],
    "Nucleotide Metabolism": [
        "nucleotide",
        "purine",
        "pyrimidine",
        "adenosine",
        "guanosine",
        "cytidine",
        "uridine",
        "thymidine",
        "inosine",
        "xanthosine",
        "deoxy",
        "dna",
        "rna",
    ],
    "Amino Acid Metabolism": [
        "amino acid",
        "aminotransferase",
        "alanine",
        "arginine",
        "aspartate",
        "asparagine",
        "cysteine",
        "glutamate",
        "glutamine",
        "glycine",
        "histidine",
        "isoleucine",
        "leucine",
        "lysine",
        "methionine",
        "ornithine",
        "phenylalanine",
        "proline",
        "serine",
        "threonine",
        "tryptophan",
        "tyrosine",
        "valine",
    ],
    "Cell Envelope Biosynthesis": [
        "cell envelope",
        "cell wall",
        "peptidoglycan",
        "lipopolysaccharide",
        "lipid a",
        "murein",
        "outer membrane",
        "porin",
        "envelope",
        "capsule",
    ],
    "Cofactor and Prosthetic Group Biosynthesis": [
        "cofactor",
        "prosthetic",
        "biotin",
        "thiamine",
        "riboflavin",
        "folate",
        "cobalamin",
        "heme",
        "porphyrin",
        "molybdopterin",
        "ubiquinone",
        "menaquinone",
        "quinone",
        "pyridox",
        "coenzyme a",
        "nad ",
        "nadp",
        "fad",
        "fmn",
    ],
    "Lipid Metabolism": [
        "lipid",
        "fatty acid",
        "acyl",
        "acp",
        "phospholipid",
        "cardiolipin",
        "glycerolipid",
        "beta oxidation",
        "sterol",
        "phosphatidyl",
        "diacylglycerol",
    ],
    "Energy Metabolism": [
        "atp synthase",
        "electron transport",
        "respiratory chain",
        "oxidative phosphorylation",
        "cytochrome",
        "nadh dehydrogenase",
        "succinate dehydrogenase",
        "terminal oxidase",
        "hydrogenase",
        "ferredoxin",
    ],
    "Carbohydrate Metabolism": [
        "carbohydrate",
        "glycolysis",
        "gluconeogenesis",
        "glucose",
        "fructose",
        "galactose",
        "sucrose",
        "maltose",
        "glycogen",
        "hexose",
        "mannose",
        "starch",
    ],
    "Alternate Carbon Metabolism": [
        "formate",
        "methanol",
        "acetate",
        "ethanol",
        "propionate",
        "butyrate",
        "glyoxylate",
        "one-carbon",
        "c1",
    ],
}

CATEGORY_PRIORITY = [
    "Photosynthesis",
    "Carotenoid Biosynthesis",
    "Aromatic Compounds Metabolism",
    "Pentose Phosphate Pathway",
    "Pyruvate Metabolism",
    "Nitrogen Metabolism",
    "Nucleotide Metabolism",
    "Amino Acid Metabolism",
    "Cell Envelope Biosynthesis",
    "Cofactor and Prosthetic Group Biosynthesis",
    "Lipid Metabolism",
    "Energy Metabolism",
    "Carbohydrate Metabolism",
    "Alternate Carbon Metabolism",
]


def parse_genes_from_gpr(gpr: str) -> List[str]:
    return [token for token in GENE_TOKEN_PATTERN.findall(gpr or "") if token.lower() not in LOGIC_WORDS]


def reaction_compartments(metabolite_ids: Iterable[str]) -> Set[str]:
    comps: Set[str] = set()
    for met_id in metabolite_ids:
        if "_" not in met_id:
            continue
        comps.add(met_id.rsplit("_", 1)[-1])
    return comps


def load_genome_annotations(gb_path: Path) -> Dict[str, Dict[str, str]]:
    annotations: Dict[str, Dict[str, str]] = {}
    for record in SeqIO.parse(str(gb_path), "genbank"):
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
            symbol = str(qualifiers.get("gene", [""])[0]).strip() if qualifiers.get("gene") else ""
            product = str(qualifiers.get("product", [""])[0]).strip() if qualifiers.get("product") else ""
            annotations[gene_id] = {"symbol": symbol, "product": product}
    return annotations


def contains_any(text: str, keywords: List[str]) -> bool:
    return any(keyword in text for keyword in keywords)


def classify_reaction(
    reaction: Dict[str, object],
    genome_annotations: Dict[str, Dict[str, str]],
) -> Tuple[str, str]:
    rid = str(reaction.get("id", "") or "")
    rid_upper = rid.upper()
    rname = str(reaction.get("name", "") or "")
    gpr = str(reaction.get("gene_reaction_rule", "") or "")
    genes = parse_genes_from_gpr(gpr)

    symbols: List[str] = []
    products: List[str] = []
    for gene_id in genes:
        ann = genome_annotations.get(gene_id)
        if ann is None:
            continue
        if ann.get("symbol"):
            symbols.append(ann["symbol"])
        if ann.get("product"):
            products.append(ann["product"])

    genome_text = " ".join([rid, rname, " ".join(symbols), " ".join(products)]).lower()
    reaction_text = f"{rid} {rname}".lower()

    metabolites = reaction.get("metabolites", {})
    if not isinstance(metabolites, dict):
        metabolites = {}
    compartments = reaction_compartments(metabolites.keys())

    if rid_upper.startswith("EX_") or rid_upper.startswith("DM_") or rid_upper.startswith("SK_"):
        return "Extracellular exchange", "model_structure"

    if len(compartments) > 1:
        return "Transport", "model_structure"

    if "transport" in rname.lower():
        return "Transport", "model_structure"

    if any(marker in rid for marker in ["TEX", "T2PP", "T3PP", "ABCPP", "PPt", "TR_"]):
        return "Transport", "model_structure"

    if products or symbols:
        for category in CATEGORY_PRIORITY:
            if contains_any(genome_text, KEYWORDS[category]):
                return category, "genome_annotation_keyword"

    for category in CATEGORY_PRIORITY:
        if contains_any(reaction_text, KEYWORDS[category]):
            return category, "reaction_keyword_fallback"

    return "Others", "unassigned"


def write_table(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "category",
                "reaction_count",
                "percent",
                "fraction",
                "total_reactions",
                "basis_breakdown",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    model_json = json.loads(FINAL_MODEL.read_text(encoding="utf-8"))
    genome_annotations = load_genome_annotations(GENOME_GB)

    category_counts: Counter[str] = Counter()
    category_basis_counts: Dict[str, Counter[str]] = {category: Counter() for category in CATEGORY_ORDER}

    reactions = model_json.get("reactions", [])
    if not isinstance(reactions, list):
        reactions = []

    for reaction in reactions:
        if not isinstance(reaction, dict):
            continue
        category, basis = classify_reaction(reaction, genome_annotations)
        category_counts[category] += 1
        category_basis_counts[category][basis] += 1

    total = sum(category_counts.values())
    rows: List[Dict[str, object]] = []
    for category in CATEGORY_ORDER:
        count = category_counts.get(category, 0)
        fraction = (count / total) if total else 0.0
        basis_breakdown = " | ".join(
            f"{basis}:{value}" for basis, value in category_basis_counts[category].most_common()
        )
        rows.append(
            {
                "category": category,
                "reaction_count": count,
                "percent": round(fraction * 100.0, 3),
                "fraction": round(fraction, 6),
                "total_reactions": total,
                "basis_breakdown": basis_breakdown,
            }
        )

    write_table(OUT_TABLE, rows)
    print(f"Exported pie table: {OUT_TABLE}")


if __name__ == "__main__":
    main()
