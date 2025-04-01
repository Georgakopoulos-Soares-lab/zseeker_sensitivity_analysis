#!/usr/bin/env python3

import csv
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def main(input_file):
    # Define the discrete ranges for GC, GT, and AC weights.
    gc_range = [4, 5, 6, 7, 8, 9]
    gt_range = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    ac_range = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]

    # Initialize counts for each parameter using the defined discrete values.
    gc_counts = {val: 0 for val in gc_range}
    gt_counts = {val: 0 for val in gt_range}
    ac_counts = {val: 0 for val in ac_range}
    
    # List to store rows where optimal_parameters == 1.
    optimal_rows = []

    # Open and read the file.
    with open(input_file, newline='') as f:
        # Use ',' as delimiter (change to '\t' if file is tab-delimited)
        reader = csv.reader(f, delimiter=',')
        headers = next(reader, None)  # Read header row.
        if not headers:
            print("The file appears to be empty or malformed.")
            return

        # Identify column indices from header.
        try:
            gc_idx = headers.index("GC_weight")
            gt_idx = headers.index("GT_weight")
            ac_idx = headers.index("AC_weight")
            opt_idx = headers.index("optimal_parameters")
        except ValueError as e:
            print("Could not find required columns in header:", e)
            return

        # Process each row.
        for row in reader:
            # Check if row has enough columns.
            if len(row) <= opt_idx:
                continue

            try:
                opt_value = float(row[opt_idx])
                print(opt_value)
            except ValueError:
                continue

            if opt_value == 1:
                # Save the optimal row.
                optimal_rows.append(row)

                # Parse GC, GT, AC weights.
                try:
                    gc_val = float(row[gc_idx])
                    gt_val = float(row[gt_idx])
                    ac_val = float(row[ac_idx])
                except ValueError:
                    continue

                # Increment the counts only if the value is in the predefined ranges.
                if gc_val in gc_counts:
                    gc_counts[gc_val] += 1
                if gt_val in gt_counts:
                    gt_counts[gt_val] += 1
                if ac_val in ac_counts:
                    ac_counts[ac_val] += 1

    # --- Print the rows where optimal_parameters == 1 ---
    print("Rows with optimal_parameters == 1:")
    print(headers)  # Print header row for context.
    for row in optimal_rows:
        print(row)

    # --- Plotting ---
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    # --- GC Weight Bar Plot ---
    gc_keys = list(gc_counts.keys())
    gc_vals = list(gc_counts.values())
    axs[0].bar(gc_keys, gc_vals, width=0.6, color='skyblue')
    axs[0].set_title("GC Weight (optimal_parameters=1)")
    axs[0].set_xlabel("GC_weight")
    axs[0].set_ylabel("Count")
    axs[0].set_xticks(gc_range)
    # Set y-axis tick interval to 1.
    axs[0].yaxis.set_major_locator(ticker.MultipleLocator(1))

    # --- GT Weight Bar Plot ---
    gt_keys = list(gt_counts.keys())
    gt_vals = list(gt_counts.values())
    axs[1].bar(gt_keys, gt_vals, width=0.06, color='lightgreen')
    axs[1].set_title("GT Weight (optimal_parameters=1)")
    axs[1].set_xlabel("GT_weight")
    axs[1].set_xticks(gt_range)
    axs[1].yaxis.set_major_locator(ticker.MultipleLocator(1))

    # --- AC Weight Bar Plot ---
    ac_keys = list(ac_counts.keys())
    ac_vals = list(ac_counts.values())
    axs[2].bar(ac_keys, ac_vals, width=0.06, color='salmon')
    axs[2].set_title("AC Weight (optimal_parameters=1)")
    axs[2].set_xlabel("AC_weight")
    axs[2].set_xticks(ac_range)
    axs[2].yaxis.set_major_locator(ticker.MultipleLocator(1))

    plt.tight_layout()
    plt.savefig('barplot.png')
    print("Bar plot saved as 'barplot.png'.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_optimal_params.py <input_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    main(input_file)
