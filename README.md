# Cis-reg: Find transcription factor binding sites in genome 

This is an ongoing project of mine that will:
* Read in a user-defined score matrix for a transcription factors binding sequence
* Read in a genome ~~in blocks~~
* Find all potential binding domains matching the score matrix
* Annotate these into a .gff file, coloured by score

Future things I'd like to add: 
* Benchmarking to inform best method to process chromosomes (i.e. as blocks or string) 
* Performing a clustering analysis to identify regions enriched for a TFs binding domain
* Run and cluster for several TFs
* Make predictions based on clustering/score on what regions in genome likely regulated by TF(s)

## Options
These are well descibed in the program - run with the `-h` flag to see the following usage message:

```{perl}
*************** Cis-reg ***************
Usage: cis-reg.v3.pl [options]
--genome = genome fasta file
--matrix = score matrx file
--outfile = specify name of output file
--help
--quiet
--debug
Nick Riddiford 2016
```

Running as-is, will take in a test genome file and score matrix, and output a `.gff3` file that can be read into IGV.
