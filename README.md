# Cis-reg: Find transcription factor binding sites in genome

This is an ongoing project of mine that will:
* Read in a user-defined score matrix for the cis-regulatory element bound by a particular transcription factor
* Read in a genome ~~in blocks~~
* Find all potential binding domains in the genome for given TF (above a score threshold)
* Annotate these into a .gff file, coloured by score

## Options
These are well descibed in the program - run with the `-h` flag to see the following usage message:

```
*************** Cis-reg ***************
Usage: cis-reg.v4.3.pl [options]
--genome = genome fasta file
--matrix = score matrx file
--score  = score cuttoff for match. 1 = 100% match - not recommended to set score < 0.75
--outfile = specify name of output file
--help
--quiet
--debug
Nick Riddiford 2016
```

Running as-is, will take in a test genome file and score matrix, and output a `.gff3` file that can be read into IGV.

# To do
- [ ] Annotate both strands of DNA
- [ ] Benchmarking to inform best method to process chromosomes (i.e. as blocks or strings)
- [ ] Clustering analysis to identify regions enriched for a TF's binding domain
- [ ] Run and cluster for several TFs
- [ ] Make predictions based on clustering/score on what regions in genome likely regulated by TF(s)
- [ ] Six1 Binding domain - TCAGGTNNC (Sato et al. 2012)
