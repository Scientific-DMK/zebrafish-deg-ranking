"""Rank a supplied differential-expression table by absolute log2 fold-change.

Designed for Kaggle Notebooks. If --input is omitted, the script searches
/kaggle/input recursively and requires exactly one Excel workbook.
"""

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd


DEFAULT_INPUT_DIRECTORY = Path("/kaggle/input")
DEFAULT_OUTPUT = Path("ranked_log2FC_genes.xlsx")


def find_single_excel_file(directory: Path) -> Path:
    files = sorted(directory.rglob("*.xlsx"))
    if not files:
        raise FileNotFoundError(f"No .xlsx files were found under {directory}.")
    if len(files) > 1:
        choices = "\n".join(f"  - {path}" for path in files)
        raise RuntimeError(
            "More than one Excel file was found. Supply the intended file with "
            f"--input.\n{choices}"
        )
    return files[0]


def rank_by_absolute_log2fc(
    dataframe: pd.DataFrame,
    logfc_column: str = "logFC",
    cutoff: float = 0.0,
) -> tuple[pd.DataFrame, int]:
    if logfc_column not in dataframe.columns:
        available = ", ".join(map(str, dataframe.columns))
        raise KeyError(
            f"Column {logfc_column!r} was not found. Available columns: {available}"
        )
    if cutoff < 0:
        raise ValueError("The absolute log2FC cutoff must be zero or greater.")

    ranked = dataframe.copy()
    numeric_logfc = pd.to_numeric(ranked[logfc_column], errors="coerce")
    dropped_rows = int(numeric_logfc.isna().sum())
    ranked[logfc_column] = numeric_logfc
    ranked = ranked.dropna(subset=[logfc_column]).copy()

    ranked["abs_log2FC"] = ranked[logfc_column].abs()
    ranked = ranked[ranked["abs_log2FC"] > cutoff].copy()
    ranked["direction"] = ranked[logfc_column].map(
        lambda value: "Upregulated" if value > 0 else "Downregulated"
    )

    ranked = ranked.sort_values(
        by="abs_log2FC", ascending=False, kind="stable"
    ).reset_index(drop=True)
    ranked.insert(0, "absolute_rank", range(1, len(ranked) + 1))
    return ranked, dropped_rows


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        help="Input .xlsx file. If omitted, exactly one file must exist in /kaggle/input.",
    )
    parser.add_argument(
        "--column",
        default="logFC",
        help="Column containing signed log2 fold-change values (default: logFC).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.0,
        help="Retain genes with absolute log2FC greater than this value (default: 0).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output workbook (default: {DEFAULT_OUTPUT}).",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    input_path = args.input or find_single_excel_file(DEFAULT_INPUT_DIRECTORY)

    dataframe = pd.read_excel(input_path)
    ranked, dropped_rows = rank_by_absolute_log2fc(
        dataframe,
        logfc_column=args.column,
        cutoff=args.cutoff,
    )

    upregulated = ranked[ranked["direction"] == "Upregulated"].copy()
    downregulated = ranked[ranked["direction"] == "Downregulated"].copy()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.output) as writer:
        ranked.to_excel(writer, sheet_name="Ranked_by_abs_log2FC", index=False)
        upregulated.to_excel(writer, sheet_name="Upregulated", index=False)
        downregulated.to_excel(writer, sheet_name="Downregulated", index=False)

    ranked.to_csv(args.output.with_name("ranked_by_abs_log2FC.csv"), index=False)
    upregulated.to_csv(args.output.with_name("upregulated_genes.csv"), index=False)
    downregulated.to_csv(
        args.output.with_name("downregulated_genes.csv"), index=False
    )

    print(f"Input: {input_path}")
    print(f"Output: {args.output}")
    print(f"Ranked genes: {len(ranked)}")
    print(f"Upregulated: {len(upregulated)}")
    print(f"Downregulated: {len(downregulated)}")
    print(f"Rows removed because {args.column!r} was missing/non-numeric: {dropped_rows}")


if __name__ == "__main__":
    main()
