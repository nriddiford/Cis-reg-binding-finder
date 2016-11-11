# Cis-reg: Find transcription factor binding sites in genome 

This is an ongoing project of mine that currently aims to:
* Read in a user-defined score matrix for a transcription factors binding sequence
* Read in a genome in blocks
* Find all potential binding domains matching the score matrix
* Annotate these into a .gff file, based on score

Future things I'd like to add include: 
* Performing a clustering analysis to identify regions enriched with binding sites
* Run and cluster for several TFs
* Make predictions based on clustering/score on what regions in genome likely regulated by TF(s)

## Options
These are well descibed in the program - run with the `-h` flag to see the following usage message:

```{perl}
*************** Cis-reg ***************
Usage: cis-reg.pl [options]
--genome = genome fasta file
--matrix = score matrx file
--block-size = size of genome chuck to process in loop. Default = 10
--outfile = specify name of output file
--help
--quiet
--debug
Nick Riddiford 2016
```

## Processing chromosomes in blocks 