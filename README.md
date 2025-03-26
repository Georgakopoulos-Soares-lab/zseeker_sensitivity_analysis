# ZSeeker sensitivity test

In our sensitivity tests, we evaluated a set of sequences where a subset was designated as non–Z-DNA forming and served as negative controls, while the remaining sequences were expected to form Z-DNA. We varied the dinucleotide weights across several parameter combinations: GC weights were tested over the range [4, 5, 6, 7, 8, 9]; GT and AC weights were examined over the range [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] (with negative values prohibited); and the AT weight was maintained at 0.5 due to its nonlinear relationship with consecutive AT transitions. By applying a threshold variation of ±2 and prohibiting negative scores, our experiments identified the optimal configuration as GC = 7, GT = 1.25, AC = 1.25, and AT = 0.5. This parameter set distinguishes between sequences that form Z-DNA and those that do not, and serves as a baseline that can be adjusted to account for environmental factors such as divalent cation concentration and pH.

1) ZSeeker was tested on a series of Z-DNA forming and non–Z-DNA forming sequences in the FASTA file sensitivity.fa.

2) Sequences marked with the prefix "non" in the FASTA file are classified as non–Z-DNA forming.

3) All possible combinations of the parameter values below were tested to determine whether ZSeeker passes the sensitivity criteria. We proved that, with a threshold variation of ±2 and by prohibiting negative scores, the optimal values are GC = 7, GT = 1.25, AC = 1.25, and AT = 0.5.

4) Users can adjust the values within the given ranges to account for environmental factors such as divalent cation concentration and pH.

## Dinucleotide weights tested
GC [4, 5, 6, 7, 8, 9] # Centered at 7 ± 3 #Default at 7 
GT [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] # Default at 1.25 
AC [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] # Default at 1.25 
AT [0.5] # Fixed due to its nonlinear relationship

5) The aggregated results are provided in a summary CSV file, which includes a column indicating whether the parameter set is optimal and additional columns for each sequence in the FASTA file showing the maximum Z-DNA score produced by that run. Individual test results for each parameter combination are available in separate directories under the sweep results folder.




