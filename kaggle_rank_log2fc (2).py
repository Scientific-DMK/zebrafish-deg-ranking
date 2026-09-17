from pathlib import Path

import pandas as pd


# Find the Excel file uploaded to Kaggle
input_folder = Path("/kaggle/input")
excel_files = list(input_folder.rglob("*.xlsx"))

if len(excel_files) == 0:
    raise FileNotFoundError("No Excel file was found in /kaggle/input")

if len(excel_files) > 1:
    print("More than one Excel file was found:")
    for file in excel_files:
        print(file)
    raise ValueError("Please remove the extra files or enter the file path manually")

file_path = excel_files[0]
df = pd.read_excel(file_path)

# Column containing the signed log2 fold-change values
logfc_column = "logFC"

if logfc_column not in df.columns:
    raise KeyError(f"The column '{logfc_column}' was not found")

# The supplied DEG lists were already filtered, so no extra cutoff was used
cutoff = 0

df[logfc_column] = pd.to_numeric(df[logfc_column], errors="coerce")
df = df.dropna(subset=[logfc_column]).copy()

# Rank all genes by the size of the expression change
df["abs_log2FC"] = df[logfc_column].abs()
ranked_genes = df[df["abs_log2FC"] > cutoff].copy()
ranked_genes = ranked_genes.sort_values("abs_log2FC", ascending=False)
ranked_genes.insert(0, "Rank", range(1, len(ranked_genes) + 1))

# Keep the original sign so upregulated and downregulated genes remain separate
ranked_genes["Direction"] = ranked_genes[logfc_column].apply(
    lambda value: "Upregulated" if value > 0 else "Downregulated"
)

positive_genes = ranked_genes[ranked_genes[logfc_column] > 0].copy()
negative_genes = ranked_genes[ranked_genes[logfc_column] < 0].copy()

output_file = "ranked_log2FC_genes.xlsx"

with pd.ExcelWriter(output_file) as writer:
    ranked_genes.to_excel(writer, sheet_name="Ranked genes", index=False)
    positive_genes.to_excel(writer, sheet_name="Upregulated", index=False)
    negative_genes.to_excel(writer, sheet_name="Downregulated", index=False)

ranked_genes.to_csv("ranked_log2FC_genes.csv", index=False)

print("Finished")
print("Input file:", file_path)
print("Total genes:", len(ranked_genes))
print("Upregulated genes:", len(positive_genes))
print("Downregulated genes:", len(negative_genes))
print("Created:", output_file)
