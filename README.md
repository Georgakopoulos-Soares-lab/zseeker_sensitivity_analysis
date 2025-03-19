# ZSeeker sensitivity test

1) ZSeeker was tested in a series of forming/non forming Z-DNA  sequences in the fasta file sensitivity.fa

2) Sequences marked with the prefix "non" in the fasta file are non Z-DNA forming

3) All possible combinations of the values below were tested to see if ZSeeker passes the tests and we proved that with a threshold of +-2 in all eligible parameters of ZSeeker and prohibiting negative scores,  the values: GC=7 , GT=1.25 , AC=1.25 and AT=0.5 are the optimal. 

4) Depending on environmental factors such as divalent cations, PH etc. the users can safely adjust the values within that +-2 range  to account for those. 

## Dinucleotide weights tested
```
GC [5, 7, 9]                # 7 ± 2 
GT [1.25, 3.25]             # 1.25 + 2 (negative values are prohibited e.g -0.75)
AC [1.25, 3.25]             # 1.25 + 2 (negative values are prohibited e.g.-0.75)
AT [0.5]                    # Not variable due to nonlinear relationship with AT_consecutive_array
```       

5) The aggregated results can be found at <b>sweep_results/parameter_sweep_summary.csv</b>  and the individual tests can be found at each directory under sweep_results , named acoordingly.