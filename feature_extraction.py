from pathlib import Path
import csv

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def read_fasta(path):
    header = None
    seq_lines = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
            else:
                seq_lines.append(line)

    sequence = "".join(seq_lines).upper()
    return header, sequence

def aa_fraction(seq, aa_set):
    if not seq:
        return 0.0
    count = sum(1 for aa in seq if aa in aa_set)
    return round(count / len(seq), 4)

def amino_acid_composition(seq):
    comp = {}
    for aa in AMINO_ACIDS:
        comp[f"frac_{aa}"] = aa_fraction(seq, {aa})
    return comp

def compute_features(seq):
    features = {
        "length_aa": len(seq),
        "acidic_fraction": aa_fraction(seq, {"D", "E"}),
        "basic_fraction": aa_fraction(seq, {"K", "R", "H"}),
        "hydrophobic_fraction": aa_fraction(seq, {"A", "V", "I", "L", "M", "F", "W", "Y"}),
        "aromatic_fraction": aa_fraction(seq, {"F", "W", "Y"}),
        "glycine_fraction": aa_fraction(seq, {"G"}),
        "proline_fraction": aa_fraction(seq, {"P"}),
    }
    features.update(amino_acid_composition(seq))
    return features

project_root = Path(__file__).resolve().parent.parent
data_dir = project_root / "data"
results_dir = project_root / "results"
metadata_file = results_dir / "protein_metadata.csv"
output_file = results_dir / "protein_feature_matrix.csv"

metadata_rows = []
with open(metadata_file, "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        metadata_rows.append(row)

feature_rows = []

for row in metadata_rows:
    fasta_path = data_dir / row["filename"]
    header, seq = read_fasta(fasta_path)

    features = compute_features(seq)

    out_row = {
        "organism": row["organism"],
        "category": row["category"],
        "environment_type": row["environment_type"],
        "filename": row["filename"],
        "header": header,
    }
    out_row.update(features)
    feature_rows.append(out_row)

fieldnames = [
    "organism",
    "category",
    "environment_type",
    "filename",
    "header",
    "length_aa",
    "acidic_fraction",
    "basic_fraction",
    "hydrophobic_fraction",
    "aromatic_fraction",
    "glycine_fraction",
    "proline_fraction",
] + [f"frac_{aa}" for aa in AMINO_ACIDS]

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(feature_rows)

print("Saved:", output_file)
for row in feature_rows:
    print(
        row["organism"],
        "| length =", row["length_aa"],
        "| acidic =", row["acidic_fraction"],
        "| basic =", row["basic_fraction"],
        "| hydrophobic =", row["hydrophobic_fraction"]
    )