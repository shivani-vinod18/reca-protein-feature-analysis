from pathlib import Path
import csv
import matplotlib.pyplot as plt
import numpy as np

CORE_FEATURES = [
    "acidic_fraction",
    "basic_fraction",
    "hydrophobic_fraction",
    "aromatic_fraction",
    "glycine_fraction",
    "proline_fraction",
]

def read_feature_matrix(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def mean(values):
    if not values:
        return 0.0
    return sum(values) / len(values)

project_root = Path(__file__).resolve().parent.parent
results_dir = project_root / "results"
input_file = results_dir / "protein_feature_matrix.csv"
summary_file = results_dir / "group_feature_summary.csv"
plot_file = results_dir / "core_feature_comparison.png"

rows = read_feature_matrix(input_file)

extreme_rows = [r for r in rows if r["category"] == "extremophile"]
control_rows = [r for r in rows if r["category"] == "control"]

summary_rows = []

for feature in CORE_FEATURES:
    extreme_vals = [float(r[feature]) for r in extreme_rows]
    control_vals = [float(r[feature]) for r in control_rows]

    extreme_mean = round(mean(extreme_vals), 4)
    control_mean = round(mean(control_vals), 4)
    difference = round(extreme_mean - control_mean, 4)

    summary_rows.append({
        "feature": feature,
        "extremophile_mean": extreme_mean,
        "control_mean": control_mean,
        "difference_extreme_minus_control": difference
    })

with open(summary_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "feature",
            "extremophile_mean",
            "control_mean",
            "difference_extreme_minus_control"
        ]
    )
    writer.writeheader()
    writer.writerows(summary_rows)

features = [row["feature"] for row in summary_rows]
extreme_means = [row["extremophile_mean"] for row in summary_rows]
control_means = [row["control_mean"] for row in summary_rows]

x = np.arange(len(features))
width = 0.35

plt.figure(figsize=(10, 6))
plt.bar(x - width / 2, extreme_means, width, label="Extremophile")
plt.bar(x + width / 2, control_means, width, label="Control")
plt.xticks(x, features, rotation=30, ha="right")
plt.ylabel("Mean fraction")
plt.title("RecA Core Sequence Features: Extremophiles vs Controls")
plt.legend()
plt.tight_layout()
plt.savefig(plot_file, dpi=300)
plt.close()

print("Saved summary:", summary_file)
print("Saved plot:", plot_file)

for row in summary_rows:
    print(
        row["feature"],
        "| extremophile =", row["extremophile_mean"],
        "| control =", row["control_mean"],
        "| diff =", row["difference_extreme_minus_control"]
    )