#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def csv_directory_to_fasta(csv_dir: str, output_fasta: str):
    csv_path = Path(csv_dir)
    fasta_path = Path(output_fasta)

    if not csv_path.exists():
        print(f"Error: Directory '{csv_dir}' does not exist.")
        return

    # Collect all CSV files, excluding any combined 'summary.csv' to avoid duplicates
    csv_files = sorted(list(csv_path.glob("*.csv")))
    csv_files = [f for f in csv_files if f.stem.lower() != "summary"]

    if not csv_files:
        print(f"No individual CSV files found in '{csv_dir}'.")
        return

    print(f"Found {len(csv_files)} CSV files. Aggregating sequences...")

    sequences_written = 0

    # Open the FASTA file for writing
    with open(fasta_path, "w") as fasta_out:
        for f in csv_files:
            if "sequence" not in str(f):
                continue
            try:
                # Read the CSV file
                df = pd.read_csv(f)

                # Verify required columns exist
                if "sequence_name" not in df.columns or "sequence" not in df.columns:
                    print(f"Warning: Skipping {f.name} (missing 'sequence_name' or 'sequence' columns)")
                    continue

                # Iterate through all rows in this CSV file (handles single or multi-row CSVs)
                for _, row in df.iterrows():
                    # Handle any unexpected NaN values gracefully
                    if pd.isna(row["sequence_name"]) or pd.isna(row["sequence"]):
                        continue

                    seq_name = str(row["sequence_name"]).strip()
                    sequence = str(row["sequence"]).strip().upper()

                    # Write in FASTA standard format
                    fasta_out.write(f">{seq_name}\n{sequence}\n")
                    sequences_written += 1

            except Exception as e:
                print(f"Error processing file {f.name}: {e}")

    print(f"Success! Collected {sequences_written} sequences into: {fasta_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collect sequences from a directory of metric CSVs into a single FASTA file."
    )
    parser.add_argument(
        "--csv-dir",
        type=str,
        default="metrics", # Defaults to your script's output directory
        help="Path to the directory containing the CSV files."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="binders_50.fasta",
        help="Path where the final FASTA file should be saved."
    )

    args = parser.parse_args()
    csv_directory_to_fasta(args.csv_dir, args.output)
