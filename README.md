# Zebrafish DEG ranking

This repository contains the Python script used for downstream processing of
differentially expressed gene datasets in the zebrafish domestication study.

The script:

- imports an Excel file containing signed log2 fold-change values;
- calculates absolute log2 fold-change;
- ranks genes from highest to lowest absolute log2 fold-change;
- retains the original expression direction; and
- exports complete, upregulated and downregulated gene lists.

## Requirements

- Python 3
- pandas
- openpyxl

## Input

The input Excel file must contain a column named `logFC`. The supplied datasets
were already filtered for differential expression before being processed by
this script.

## Output

The script produces:

- `ranked_log2FC_genes.xlsx`
- `ranked_log2FC_genes.csv`

The Excel workbook contains separate sheets for all ranked genes, upregulated
genes and downregulated genes.
