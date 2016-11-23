#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature qw/ say /; 
use Getopt::Long qw/ GetOptions /;

use Bio::SeqIO;


my $genome_file = 'in.fa';
my $matrix_file = 'matrix_for_Dros_so.txt';
my $out_file = 'binding_domains';
my $threshold = 0.8;
my $debug;
my $quiet;
my $help;

# Should add score threshold option
GetOptions( 'genome=s'     =>    \$genome_file,
		   	'matrix=s'     =>    \$matrix_file,
			'score=i'	   =>	 \$threshold,
		   	'outfile=s'    =>	 \$out_file,
		   	'help'         =>    \$help,
		   	'quiet'        =>    \$quiet,
           	'debug'        =>    \$debug
) or die usage();

if ($help)  { exit usage() } 
if ($quiet) { say "Running in quiet mode" }
if ($debug) { say "Running in debug mode" }

my $start_run = time();

# need to make sure this handles errors properly
say "Reading genome file... '$genome_file'";
my $seqio = Bio::SeqIO->new(-file => $genome_file, '-format' => 'Fasta');

say "Reading matrix file... '$matrix_file'";
open my $matrix, '<', $matrix_file or die $!;
read_matrix($matrix);

say "Writing output to '$out_file'";

my (@consensus, %score_lookup);

my $domain_length = scalar(@consensus);

say "Looking for binding domain $domain_length long" if $debug;
say "Using score cutoff of >= $threshold";

my $factor = $domain_length * 4;
my %candidates;
my %kmer_freq;

open my $out, '>', "$out_file\.gff3" or die $!;

print $out "##gff-version 3\n";
print $out "#track name=\"$out_file\" color=#FFBB33 gffTags=on\n";

# Add more colours for longer elements
my $dark_red = '#900C3F';
my $red = '#FF5733';
my $orange = '#F8C211';

# How does this read-in stragegy cover line-wrapped instances? 
# Does it catch elements spanning two lines? 
while(my $seq = $seqio->next_seq) {
	my $nucs = $seq->seq;
	my $chr = $seq->id;
	my $chr_length =$seq->length; # 1 based length of chromosome
	say "Processing '$chr'...";
	my $dom_count = 0;

	for my $i ( 0 .. $chr_length - $domain_length ){
		my $kmer = substr($nucs, $i, $domain_length);
		 
		# Make 1-based for gff3
		my $start = $i + 1;
		my $stop = $i + $domain_length;
 	    
		# $kmer_freq{$kmer}++; # populate hash with frequency of k-mer
		 		 
		my @nucs = split '', $kmer;
		my $nuc_score = 0;
		 
		# iterate over each nucleotide and sum score from lookup kmer
		for ( 0 .. $#nucs ){
			if ( exists $score_lookup{$_}{$nucs[$_]} ){
		 		$nuc_score += $score_lookup{$_}{$nucs[$_]};
		  	}
		}
		
		# annotate hits with high score
		# Need to include some with less than perfect scores...
		my $score = $nuc_score/$factor;
		
		if ( $score >= $threshold ){ 
			my ($sc) = sprintf('%.0f', ($score * 100) );
		
			print "$kmer = $nuc_score  HIT!  $sc\% score against consensus\n" if $debug;
			$dom_count++;

			my $colour = $orange;
		
			if ( $score == 1 ) { $colour = $dark_red }
			if ( $score >= 0.8 and $score < 1 ) { $colour = $red }
			
			my $gff_line = join("\t", "$chr", ".", "region", "$start", "$stop", ".", "+", ".", "Name=$sc;colour=$colour");
			print $out "$gff_line\n";
			
			$candidates{$kmer}{$nuc_score}++;
		}
					
 	}
	say "Found $dom_count SO binding domains on '$chr'";
}

sub read_matrix {
	while(<$matrix>){
		my @parts = split;
		
		if ($. == 1){
			push @consensus, @parts;
		}
		else {
			for (my $cell = 0; $cell < $#parts; $cell++) {
				$score_lookup{$cell}{$parts[0]} = $parts[$cell+1];
			}
		}
	}
}
	

sub usage {
    say "*************** Cis-reg ***************";
    say "Usage: $0 [options]";
	say "--genome = genome fasta file";
	say "--matrix = score matrx file";
	say "--score  = score cuttoff for match. 1 = 100% match - not recommended to set score < 0.75";
	say "--outfile = specify name of output file";
	say "--help";
	say "--quiet";
	say "--debug";
	say "Nick Riddiford 2016";
}