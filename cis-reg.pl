#!/usr/bin/perl
use warnings;
use strict;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature qw(say);

my $start_run = time();

open my $matrix, '<', 'matrix.txt' or die $!;

my (@consensus, %score_lookup);

read_matrix($matrix);

my $infile = $ARGV[0];

open my $genome, '<', $infile or die $!;


# Six1 consensus = 'ATCCTGA';

my $binding_dom = scalar(@consensus);

use constant BLOCK_SIZE => 10;

# set params for while loop
my $block;
my $length = 0;

my $block_count = 0;
my %candidates;

open my $out, '>', 'out.txt' or die $!;

# Read into $block 1MB chunks of data from $in
# read FILEHANDLE, SCALAR, LENGTH, OFFSET
while ( my $nucs = sysread $genome, $block, BLOCK_SIZE, $length ) {

	$length += $nucs;
	
	$block_count++;
		
	# print "\nBlock: $block_count\nBlock sequence: $block\n";

	for my $offset ( 0 .. $length - $binding_dom ) {
		my $kmer = substr $block, $offset, $binding_dom;
				
		my @nucs = split '', $kmer;
		my $nuc_score = 0;
				
		# iterate over each nucleotide and sum score from lookup kmer
		for ( 0 .. $#nucs ){
			if ( exists $score_lookup{$_}{$nucs[$_]} ){
				$nuc_score += $score_lookup{$_}{$nucs[$_]};	
			}
		}	
		
		# annotate hits with score >= 10
		if ( $nuc_score >= 12 ){
			my ($sc) = (($nuc_score/16));
			print "$kmer = $nuc_score  HIT!  $sc score against consensus\n";
			print $out "||Six1||$kmer=$sc|";
			$candidates{$kmer}{$nuc_score}++;
		}

		else {
			# print "$kmer = $nuc_score\n";
			# print $out "$block\n";
		}
		
    }
	my $re_print = $block;
	
	$re_print = substr $block, ($binding_dom-1) unless $block_count == 1; # this ensures that the blocks are non-overlapping
	
	print $out "$re_print";
	
    $block = substr $block, - ($binding_dom-1);
	$length = length $block;
}
my ($seq) = join('', @consensus);
say "Top matches to $seq";

for my $match (keys %candidates ){
	for my $score ( sort { $candidates{$b} <=> $candidates{$a} } keys $candidates{$match} ){
		print "$match, $score = $candidates{$match}{$score}\n";
	}
}

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Script took $run_time seconds to complete\n";

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