#!/usr/bin/env python3
import subprocess
import itertools
import os
import pandas as pd
from pathlib import Path
import argparse
import re

def get_fasta_records(fasta_path):
    """
    Parse the FASTA file and return a dictionary mapping record IDs to sequences.
    Assumes each record's header is of the form '>recordID ...' and that the sequence lines follow.
    """
    records = {}
    with open(fasta_path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records[header] = ''.join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records[header] = ''.join(seq_lines)
    return records

def run_zseeker(fasta, gc_weight, gt_weight, ac_weight, at_weight, output_dir):
    """
    Runs the ZSeeker command-line tool with the given parameters.
    Creates a unique output directory for each parameter set.
    Returns the path to the resulting CSV file.
    """
    unique_output = Path(output_dir) / f"GC{gc_weight}_GT{gt_weight}_AC{ac_weight}_AT{at_weight}"
    unique_output.mkdir(parents=True, exist_ok=True)

    cmd = [
        "ZSeeker",
        "--fasta", str(fasta),
        "--GC_weight", str(gc_weight),
        "--GT_weight", str(gt_weight),
        "--AC_weight", str(ac_weight),
        "--AT_weight", str(at_weight),
        "--total_sequence_scoring",
        "--output_dir", str(unique_output)
    ]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running ZSeeker with GC:{gc_weight}, GT:{gt_weight}, AC:{ac_weight}, AT:{at_weight}")
        print(result.stderr)

    # The generated CSV is typically named <fasta.stem>_zdna_score.csv
    return unique_output / f"{Path(fasta).stem}_zdna_score.csv"

def evaluate_output(csv_file):
    """
    Reads the Z-DNA output CSV (which now has one row per sequence) and checks:
      - If the Chromosome (after stripping spaces and lowering) contains the postfix "_non",
        it is considered NOT Z-DNA forming and must have a Z-DNA Score < 50.
      - Otherwise, it is considered Z-DNA forming and must have a Z-DNA Score > 50.
    Returns:
      - run_df: a DataFrame of the CSV contents.
      - optimal: 1 if the parameter set meets the conditions for every sequence, 0 otherwise.
    """
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return pd.DataFrame(), 0

    if df.empty:
        return df, 0

    optimal = 1
    for idx, row in df.iterrows():
        chrom = str(row['Chromosome']).strip()
        try:
            score = float(str(row['Z-DNA Score']).strip())
        except Exception:
            score = None

        # Check for the _non postfix (case-insensitive)
        if chrom.lower().endswith("_non"):
            # For control sequences, score must be strictly less than 50.
            if score is None or score >= 50:
                optimal = 0
                break
        else:
            # For Z-DNA forming sequences, score must be strictly greater than 50.
            if score is None or score <= 50:
                optimal = 0
                break

    return df, optimal

def get_max_scores_per_sequence(run_df, seq_ids):
    """
    Given the run_df with columns [Chromosome, Start, End, Z-DNA Score, Sequence],
    extract the maximum Z-DNA Score for each sequence ID in seq_ids.
    Returns a dict {seqID: maxScore or '' if none}.
    
    - Since there is one row per sequence, the value is just the reported score.
    """
    run_df["Z-DNA Score"] = pd.to_numeric(run_df["Z-DNA Score"], errors="coerce")
    max_scores = {}
    for seq_id in seq_ids:
        subset = run_df[run_df["Chromosome"] == seq_id]
        if not subset.empty:
            max_score = subset["Z-DNA Score"].max()
            max_scores[seq_id] = max_score if not pd.isna(max_score) else ''
        else:
            max_scores[seq_id] = ''
    return max_scores

def main():
    parser = argparse.ArgumentParser(
        description="Parameter sweep for ZSeeker with expanded parameter ranges. "
                    "Generates a CSV summarizing each run, indicating whether it's optimal, "
                    "and storing the maximum Z-DNA score for each sequence in the FASTA."
    )
    parser.add_argument("--fasta", type=str, required=True, help="Path to the FASTA file (e.g., sensitivity.fa)")
    parser.add_argument("--output_dir", type=str, default="sweep_results", help="Directory for storing run outputs and summary")
    args = parser.parse_args()

    fasta = Path(args.fasta).expanduser().resolve()
    if not fasta.is_file():
        print("Error: FASTA file does not exist.")
        return

    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Parse the FASTA file to get a mapping of sequence IDs to their sequences.
    records = get_fasta_records(fasta)
    seq_ids = list(records.keys())

    # Define the updated parameter ranges.
    gc_range = [4, 5, 6, 7, 8, 9]  # 7 ± 2 => expanded
    gt_range = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    ac_range = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    at_range = [0.5]

    summary_rows = []

    for gc, gt, ac, at in itertools.product(gc_range, gt_range, ac_range, at_range):
        print(f"Running parameters: GC={gc}, GT={gt}, AC={ac}, AT={at}")
        output_csv = run_zseeker(fasta, gc, gt, ac, at, base_output_dir)
        run_df, optimal = evaluate_output(output_csv)

        row_data = {
            "GC_weight": gc,
            "GT_weight": gt,
            "AC_weight": ac,
            "AT_weight": at,
            "optimal_parameters": optimal
        }

        if not run_df.empty:
            scores_map = get_max_scores_per_sequence(run_df, seq_ids)
        else:
            scores_map = {sid: '' for sid in seq_ids}

        row_data.update(scores_map)
        summary_rows.append(row_data)

    # Build DataFrame with fixed columns and one column per sequence (using original IDs).
    columns = ["GC_weight", "GT_weight", "AC_weight", "AT_weight", "optimal_parameters"] + seq_ids
    summary_df = pd.DataFrame(summary_rows, columns=columns)

    # Rename sequence columns to include the sequence string in the header.
    rename_dict = {seq_id: f"{seq_id}( {records[seq_id]} )" for seq_id in seq_ids}
    summary_df.rename(columns=rename_dict, inplace=True)

    summary_csv = base_output_dir / "parameter_sweep_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Parameter sweep completed. Summary written to {summary_csv}")

if __name__ == "__main__":
    main()
