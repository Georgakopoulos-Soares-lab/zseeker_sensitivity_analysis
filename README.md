# ZSeeker sensitivity test

In our sensitivity tests, we evaluated a set of sequences where a subset was designated as non–Z-DNA forming and served as negative controls, while the remaining sequences were expected to form Z-DNA. We varied the dinucleotide weights across several parameter combinations: GC weights were tested over the range [4, 5, 6, 7, 8, 9]; GT and AC weights were examined over the range [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] (with negative values prohibited); and the AT weight was maintained at 0.5 due to its nonlinear relationship with consecutive AT transitions.
1) ZSeeker was tested on a series of Z-DNA forming and non–Z-DNA forming sequences in the FASTA file sensitivity.fa.

2) Sequences marked with the prefix "non" in the FASTA file are classified as non–Z-DNA forming.

3) All possible combinations of the parameter values below were tested to determine whether ZSeeker passes the sensitivity criteria. We proved that,  the optimal values are GC = 7, GT = 1.25 or 1, AC = 1.25 or 1, and AT = 0.5.

4) Users can adjust the values within the given ranges to account for environmental factors such as divalent cation concentration and pH.

## Dinucleotide weights tested
GC [4, 5, 6, 7, 8, 9] # Centered at 7 ± 3 #Default at 7 
GT [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] # Default at 1.25 
AC [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] # Default at 1.25 
AT [0.5] # Fixed due to its nonlinear relationship

5) The aggregated results are provided in a summary CSV file, which includes a column indicating whether the parameter set is optimal and additional columns for each sequence in the FASTA file showing the maximum Z-DNA score produced by that run. Individual test results for each parameter combination are available in separate directories under the sweep results folder.

6) To re-run the sensitivity analysis users have to issue the following python command:
```python
python3 sensitivity.py --fasta sensitivity.fa
```

The aggregated results are stored in:
```
 sweep_results/parameter_sweep_summary.csv
```

# Sequences used for testing

All data the were used in the sensitivity analysis were derived from wet-lab experiments and were enriched with sequences derived from the original paper of ZHunt:  
[Schroth GP, Chou PJ, Ho PS. Mapping Z-DNA in the human genome. Computer-aided mapping reveals a nonrandom distribution of potential Z-DNA-forming sequences in human genes. J Biol Chem 1992;](https://pubmed.ncbi.nlm.nih.gov/1601856/)
Sequences marked with the postfix "_non" are marked as non Z-DNA forming otherwise they are marked as forming.
To be marked with the "optimal_parameters" property a test has to predict correctly all sequences in the test.


# Analyzing Results

To analyze the results produced by the sensitivity analysis script we used a python script that plots the times each parameter value is present in in a parameter set that is optimal and by optimal we mean that it predicts all  sequences that are forming as Z-DNA forming (threshold > 50) and all sequences that are not forming as not Z-DNA forming (threshold<50).  

The analysis clearly demonstrates that a GC weight of 7 is optimal and GT and AC weights of 1 and 1.25 are also optimal by having the same count.Despite the added complexity of two extra decimal points we chose 1.25 because that aligns better with the results derived from Z-DNA formation experiments

[![barplot-2.png](https://i.postimg.cc/VNRGD2bL/barplot-2.png)](https://postimg.cc/XZrkvH8T)

The results of the analysis can be reproduced by running
```python
python barplot.py ./sweep_results/parameter_sweep_summary.csv
```




