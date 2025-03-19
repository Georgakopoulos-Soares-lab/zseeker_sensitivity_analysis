#!/usr/bin/env python3
import subprocess
import itertools
import os
import pandas as pd
from pathlib import Path
import argparse

def run_zseeker(fasta, gc_weight, gt_weight, ac_weight, at_weight, output_dir):
    """
    Runs the ZSeeker command-line tool with the given parameters.
    The output_dir is set uniquely for each run.
    """
    # Construct a unique output directory for this parameter set.
    unique_output = Path(output_dir) / f"GC{gc_weight}_GT{gt_weight}_AC{ac_weight}_AT{at_weight}"
    unique_output.mkdir(parents=True, exist_ok=True)
    # Build the command. (Other parameters can be added as needed.)
    cmd = [
        "ZSeeker",
        "--fasta", str(fasta),
        "--GC_weight", str(gc_weight),
        "--GT_weight", str(gt_weight),
        "--AC_weight", str(ac_weight),
        "--AT_weight", str(at_weight),
        "--output_dir", str(unique_output)
    ]
    # Run the command and wait for completion.
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running ZSeeker with GC:{gc_weight}, GT:{gt_weight}, AC:{ac_weight}, AT:{at_weight}")
        print(result.stderr)
    return unique_output / f"{Path(fasta).stem}_zdna_score.csv"

def evaluate_output(csv_file):
    """
    Reads the Z-DNA output CSV and checks the following conditions:
      - For every sequence whose header contains "non", there is no Z-DNA score.
      - For every sequence that does NOT contain "non", the Z-DNA score is > 50.
    Returns:
      - A DataFrame with one row per subarray from the CSV (if any)
      - A flag (1 if the parameter set is optimal, 0 otherwise)
    """
    # Read CSV. Expect columns: Chromosome, Start, End, Z-DNA Score, Sequence
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return pd.DataFrame(), 0

    # Standardize: if there is no Z-DNA score (NaN or blank) we treat it as missing.
    # For optimal parameters:
    #   All sequences with "non" in Chromosome should have no valid score (NaN)
    #   All others should have a numeric score > 50.
    optimal = 1
    # For each unique sequence (using the Chromosome field as id)
    # Here we assume that each row represents one candidate subarray.
    # For checking conditions, we group by Chromosome.
    groups = df.groupby("Chromosome")
    for name, group in groups:
        if "non" in name.lower():
            # For non sequences, no subarray should have a valid score.
            # If any row has a non-null score, mark as non-optimal.
            if group["Z-DNA Score"].dropna().shape[0] > 0:
                optimal = 0
                break
        else:
            # For others, every row should have a numeric score > 50.
            scores = pd.to_numeric(group["Z-DNA Score"], errors="coerce")
            if scores.dropna().empty or (scores.dropna() <= 50).any():
                optimal = 0
                break
    return df, optimal

def main():
    parser = argparse.ArgumentParser(
        description="Parameter sweep for ZSeeker. Runs the ZSeeker CLI on a FASTA file using various parameter combinations, "
                    "evaluates the output, and produces a summary CSV with the parameters and Z-DNA scores."
    )
    parser.add_argument("--fasta", type=str, required=True, help="Path to the FASTA file (e.g., sensitivity.fa)")
    parser.add_argument("--output_dir", type=str, default="sweep_results", help="Directory for storing run outputs and summary")
    args = parser.parse_args()

    fasta = Path(args.fasta).expanduser().resolve()
    if not fasta.is_file():
        print(f"Error: FASTA file {fasta} does not exist.")
        return

    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Define parameter ranges:
    gc_range = [5, 7, 9]                # 7 ± 2 
    gt_range = [1.25, 3.25]         # 1.25 + 2 (negative values are prohibited e.g -0.75)
    ac_range = [1.25, 3.25]         # 1.25 + 2 (negative values are prohibited e.g.-0.75)
    at_range = [0.5]            # 0.5
    # List to collect summary rows
    summary_rows = []

    # Loop over all combinations (81 total)
    for gc, gt, ac, at in itertools.product(gc_range, gt_range, ac_range, at_range):
        print(f"Running parameters: GC={gc}, GT={gt}, AC={ac}, AT={at}")
        output_csv = run_zseeker(fasta, gc, gt, ac, at, base_output_dir)
        # Evaluate output from this run:
        run_df, optimal = evaluate_output(output_csv)
        # For each subarray reported in the run_df, record a row in the summary.
        # If no subarray is reported, we still record one row with empty score.
        if run_df.empty:
            # Use filename (or parameter set id) as sequence identifier placeholder.
            summary_rows.append({
                "sequence": "",
                "score": "",
                "GC_weight": gc,
                "GT_weight": gt,
                "AC_weight": ac,
                "AT_weight": at,
                "optimal_parameters": optimal
            })
        else:
            for _, row in run_df.iterrows():
                summary_rows.append({
                    "sequence": row.get("Sequence", ""),
                    "score": row.get("Z-DNA Score", ""),
                    "GC_weight": gc,
                    "GT_weight": gt,
                    "AC_weight": ac,
                    "AT_weight": at,
                    "optimal_parameters": optimal
                })

    # Create summary dataframe
    summary_df = pd.DataFrame(summary_rows, columns=["sequence", "score", "GC_weight", "GT_weight", "AC_weight", "AT_weight", "optimal_parameters"])
    summary_csv = base_output_dir / "parameter_sweep_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Parameter sweep completed. Summary written to {summary_csv}")

if __name__ == "__main__":
    main()
