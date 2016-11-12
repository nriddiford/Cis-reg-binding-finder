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
my $block_size = 10;
my $out_file = 'SO_doms.gff3';
my $debug;
my $quiet;
my $help;

GetOptions( 'genome=s'      =>   \$genome_file,
		   	'matrix=s'     =>    \$matrix_file,
		   	'block-size=i' =>    \$block_size,
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

my (@consensus, %score_lookup);

my $binding_dom = scalar(@consensus);

my %candidates;
my %data;

open my $out, '>', $out_file or die $!;

print $out "##gff-version 3\n";
print $out "#track name=\"SO binding domains\" color=#FFBB33 gffTags=on\n";

while(my $seq = $seqio->next_seq) {
	my $nucs = $seq->seq;
	my $chr = $seq->id;
	my $chr_length =$seq->length; # 1 based length of chromosome

	for my $i ( 0 .. $chr_length - $binding_dom ) {
 	     my $kmer = substr($nucs, $i, $binding_dom);
		 
		 my $start = $i + 1;
		 my $stop = $i + 6;
 	    
		$data{$kmer}++; # populate hash with frequency of k-mer
		 		 
		 my @nucs = split '', $kmer;
		 my $nuc_score = 0;
		 
		 # iterate over each nucleotide and sum score from lookup kmer
		 	for ( 0 .. $#nucs ){
		 		if ( exists $score_lookup{$_}{$nucs[$_]} ){
		 			$nuc_score += $score_lookup{$_}{$nucs[$_]};
		  		}
		  	}
			
			# annotate hits with score >= 10
			if ( $nuc_score >= 20 ){ # need to adjust this section so that the score is automatically calculated depending on the number of nucs in consensus
				my ($sc) = (($nuc_score/24));

				 # "$kmer = $nuc_score  HIT!  $sc score against consensus\n" if $debug;
				my $gff_line = join("\t", "$chr", ".", "region", "$start", "$stop", ".", "+", ".", "Name=$nuc_score;color=#FFBB33");
				print $out "$gff_line\n";
				
				# X	.	region	3134870	3172221	.	+	.	Name=Notch;color=#FFBB33
				
				# join("\t", $chr \. region);#$i\t$stop\t\.\t\+\t\.\tName=SO_\$nuc_score;color=#FFBB33;
				$candidates{$kmer}{$nuc_score}++;
			}

			else {
				# print "$kmer = $nuc_score\n";
				# print $out "$block\n";
			}
					
 	}
	

}





 






# Read into $block 1MB chunks of data from $in
# read FILEHANDLE, SCALAR, LENGTH, OFFSET
# while ( my $nucs = sysread $genome, $block, $block_size, $length ) {
#
# 	if ($debug){
# 		print "nucs = $nucs\n";
# 		print "block = $block\n";
# 		print "block size = $block_size\n";
# 		print "length = $length\n";
# 	}
# 	$length += $nucs;
#
# 	$block_count++;
#
# 	# print "\nBlock: $block_count\nBlock sequence: $block\n";
#
# 	for my $offset ( 0 .. $length - $binding_dom ) {
# 		my $kmer = substr $block, $offset, $binding_dom;
#
# 		my @nucs = split '', $kmer;
# 		my $nuc_score = 0;
#
# 		# iterate over each nucleotide and sum score from lookup kmer
# 		for ( 0 .. $#nucs ){
# 			if ( exists $score_lookup{$_}{$nucs[$_]} ){
# 				$nuc_score += $score_lookup{$_}{$nucs[$_]};
# 			}
# 		}
#
# 		# annotate hits with score >= 10
# 		if ( $nuc_score >= 22 ){ # need to adjust this section so that the score is automaticall calculated depending on the number of nucs in consensus
# 			my ($sc) = (($nuc_score/24));
#
# 			print "$kmer = $nuc_score  HIT!  $sc score against consensus\n" unless $quiet;
# 			print $out "||Six1||$kmer=$sc|" unless $quiet;
#
# 			$candidates{$kmer}{$nuc_score}++;
# 		}
#
# 		else {
# 			# print "$kmer = $nuc_score\n";
# 			# print $out "$block\n";
# 		}
#
#     }
# 	my $re_print = $block;
#
# 	$re_print = substr $block, ($binding_dom-1) unless $block_count == 1; # this ensures that the blocks are non-overlapping
#
# 	print $out "$re_print" unless $quiet;
#
#     $block = substr $block, - ($binding_dom-1);
# 	$length = length $block;
# }
# my ($seq) = join('', @consensus);
# say "Top matches to $seq" unless $quiet;
#
# for my $match (keys %candidates ){
# 	for my $score ( sort { $candidates{$b} <=> $candidates{$a} } keys $candidates{$match} ){
# 		print "$match, $score = $candidates{$match}{$score}\n" unless $quiet;
# 	}
# }
#
# my $end_run = time();
# my $run_time = $end_run - $start_run;
# print "Script ran in $run_time seconds\n";

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
	say "--block-size = size of genome chuck to process in loop. Default = $block_size";
	say "--outfile = specify name of output file";
	say "--help";
	say "--quiet";
	say "--debug";
	say "Nick Riddiford 2016";
}