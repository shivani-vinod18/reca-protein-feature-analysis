from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

def read_csv(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

project_root = Path(__file__).resolve().parent.parent
results_dir = project_root / "results"
input_file = results_dir / "protein_feature_matrix.csv"
coords_file = results_dir / "pca_coordinates.csv"
plot_file = results_dir / "pca_plot.png"
summary_file = results_dir / "pca_summary.txt"

rows = read_csv(input_file)

# Use numeric features only
feature_names = [
    "length_aa",
    "acidic_fraction",
    "basic_fraction",
    "hydrophobic_fraction",
    "aromatic_fraction",
    "glycine_fraction",
    "proline_fraction",
] + [f"frac_{aa}" for aa in list("ACDEFGHIKLMNPQRSTVWY")]

organisms = [row["organism"] for row in rows]
categories = [row["category"] for row in rows]

X = np.array([[float(row[f]) for f in feature_names] for row in rows], dtype=float)

# Standardize features
means = X.mean(axis=0)
stds = X.std(axis=0)
stds[stds == 0] = 1.0
X_scaled = (X - means) / stds

# PCA using numpy SVD
U, S, VT = np.linalg.svd(X_scaled, full_matrices=False)
components = VT[:2]
scores = X_scaled @ components.T

# Explained variance ratio
eigenvalues = (S ** 2) / (X_scaled.shape[0] - 1)
explained_variance_ratio = eigenvalues / eigenvalues.sum()
pc1_var = explained_variance_ratio[0]
pc2_var = explained_variance_ratio[1]

# Save coordinates
with open(coords_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["organism", "category", "PC1", "PC2"])
    for org, cat, coord in zip(organisms, categories, scores):
        writer.writerow([org, cat, round(coord[0], 4), round(coord[1], 4)])

# Plot
plt.figure(figsize=(9, 7))

for cat in sorted(set(categories)):
    idx = [i for i, c in enumerate(categories) if c == cat]
    plt.scatter(
        scores[idx, 0],
        scores[idx, 1],
        label=cat,
        s=70
    )
    for i in idx:
        plt.text(scores[i, 0], scores[i, 1], organisms[i], fontsize=8)

plt.xlabel(f"PC1 ({pc1_var:.2%} variance)")
plt.ylabel(f"PC2 ({pc2_var:.2%} variance)")
plt.title("PCA of RecA Protein Sequence Features")
plt.legend()
plt.tight_layout()
plt.savefig(plot_file, dpi=300)
plt.close()

# Save summary
with open(summary_file, "w") as f:
    f.write("PCA Summary\n")
    f.write(f"Number of proteins: {len(rows)}\n")
    f.write(f"Number of features: {len(feature_names)}\n")
    f.write(f"PC1 explained variance: {pc1_var:.4f}\n")
    f.write(f"PC2 explained variance: {pc2_var:.4f}\n")
    f.write(f"Combined PC1+PC2 variance: {(pc1_var + pc2_var):.4f}\n")

print("Saved:", coords_file)
print("Saved:", plot_file)
print("Saved:", summary_file)
print(f"PC1 explained variance: {pc1_var:.4f}")
print(f"PC2 explained variance: {pc2_var:.4f}")
print(f"Combined PC1+PC2 variance: {(pc1_var + pc2_var):.4f}")