#!/usr/bin/env python3
import subprocess
import itertools
import os
import pandas as pd
from pathlib import Path
import argparse
import re

def get_fasta_ids(fasta_path):
    """
    Quickly parse the FASTA file to retrieve a list of record IDs.
    Assumes each record's first line is of the form '>recordID ...'
    Returns a list of record IDs (e.g. ['Z21', 'Z22', 'Z26', ...]).
    """
    ids = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                # For example, '>Z21 something' -> 'Z21'
                header = line[1:].split()[0]
                ids.append(header)
    return ids

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
    Reads the Z-DNA output CSV and checks:
      1) For every sequence whose Chromosome contains 'non', there is no valid Z-DNA score.
      2) For every sequence that does NOT contain 'non', the Z-DNA score is > 50 (for all subarrays).
    Returns:
      - run_df: a DataFrame of the subarrays from the CSV (may be empty if no subarrays found)
      - optimal: 1 if the parameter set meets the above conditions, 0 otherwise
    """
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return pd.DataFrame(), 0

    # If there's no data, it can't pass the test if any "non" sequences appear,
    # but let's do a straightforward check.
    if df.empty:
        # If the CSV is empty, check if that ironically might pass or fail:
        # - All 'non' would have no subarrays => that's okay,
        # - But all normal sequences also need a score > 50, which isn't satisfied.
        # We'll remain consistent with the previous logic:
        return df, 0

    optimal = 1
    groups = df.groupby("Chromosome")
    for name, group in groups:
        if "non" in name.lower():
            # Should have no valid Z-DNA subarrays
            if group["Z-DNA Score"].dropna().shape[0] > 0:
                optimal = 0
                break
        else:
            # All subarrays should have score > 50
            scores = pd.to_numeric(group["Z-DNA Score"], errors="coerce")
            if scores.dropna().empty or (scores.dropna() <= 50).any():
                optimal = 0
                break

    return df, optimal

def get_max_scores_per_sequence(run_df, seq_ids):
    """
    Given the run_df with columns [Chromosome, Start, End, Z-DNA Score, Sequence],
    extract the maximum Z-DNA Score for each sequence ID in seq_ids.
    Returns a dict {seqID: maxScore or '' if none}.

    - If multiple subarrays are reported for a sequence ID, we take the max score.
    - If no subarray is reported, return '' (empty) for that seq ID.
    """
    # Convert 'Z-DNA Score' to numeric
    run_df["Z-DNA Score"] = pd.to_numeric(run_df["Z-DNA Score"], errors="coerce")

    # We'll store the max score for each Chromosome. 
    # If a Chromosome is repeated, we take the max among its subarrays.
    # Then we match Chromosome back to seq_ids if there's a match.
    max_scores = {}
    grouped = run_df.groupby("Chromosome")["Z-DNA Score"].max()
    # grouped is a Series: index=Chromosome, value=maxScore

    # Build a mapping Chromosome->maxScore
    chrom_to_max = grouped.to_dict()

    # For each seq ID from the FASTA, see if it appears in chrom_to_max
    out = {}
    for seq_id in seq_ids:
        if seq_id in chrom_to_max and not pd.isna(chrom_to_max[seq_id]):
            out[seq_id] = chrom_to_max[seq_id]
        else:
            out[seq_id] = ''
    return out

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
        print(f"Error: FASTA file does not exist.")
        return

    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Collect the list of sequence IDs from the FASTA. We'll use these as column names.
    seq_ids = get_fasta_ids(fasta)

    # Define the updated parameter ranges
    gc_range = [4, 5, 6, 7, 8, 9]  # 7 ± 2 => expanded
    gt_range = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    ac_range = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    at_range = [0.5]

    # Prepare to collect summary
    summary_rows = []

    # Iterate over all combinations
    for gc, gt, ac, at in itertools.product(gc_range, gt_range, ac_range, at_range):
        print(f"Running parameters: GC={gc}, GT={gt}, AC={ac}, AT={at}")
        output_csv = run_zseeker(fasta, gc, gt, ac, at, base_output_dir)
        run_df, optimal = evaluate_output(output_csv)

        # For each parameter set, we produce exactly one summary row.
        # That row has columns for GC_weight, GT_weight, AC_weight, AT_weight, optimal, 
        # plus a column for each sequence in seq_ids with that run's maximum Z-DNA score.
        row_data = {
            "GC_weight": gc,
            "GT_weight": gt,
            "AC_weight": ac,
            "AT_weight": at,
            "optimal_parameters": optimal
        }

        # gather max score per sequence
        if not run_df.empty:
            max_scores_map = get_max_scores_per_sequence(run_df, seq_ids)
        else:
            max_scores_map = {sid: '' for sid in seq_ids}

        # Merge these into row_data
        row_data.update(max_scores_map)

        summary_rows.append(row_data)

    # Now we form a DataFrame. The columns are 
    # ["GC_weight","GT_weight","AC_weight","AT_weight","optimal_parameters"] + seq_ids
    columns = ["GC_weight","GT_weight","AC_weight","AT_weight","optimal_parameters"] + seq_ids
    summary_df = pd.DataFrame(summary_rows, columns=columns)

    # Write out final summary
    summary_csv = base_output_dir / "parameter_sweep_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Parameter sweep completed. Summary written to {summary_csv}")

if __name__ == "__main__":
    main()
