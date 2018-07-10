#consider all inter/intra chromosomal contact
#1 remove intra_fragment
#2 remove low mapq 
#3 get fragment file 

use File::Basename;
use POSIX;
use List::Util qw[min max];
use Getopt::Std;
use vars qw/ $opt_s $opt_q /;

# Check arguments
getopts('s:q');

my $site_file = "";
my $mapq_threshold = 1;


if ($opt_s) {
  $site_file = $opt_s;
}

if ($opt_q) {
  $mapq_threshold = $opt_q;
}

# Global variables
my %chromosomes;

if (index($site_file, "none") != -1) {
   #no restriction enzyme
   die "no retriction enzyme site file!\n";
}
else {
  # read in restriction site file and store as multidimensional array
  open FILE, $site_file or die $!;
  while (<FILE>) {
    my @locs = split;
    my $key = shift(@locs);
    my $ref = \@locs;
    $chromosomes{$key} = $ref;
  }
  close(FILE);
}

# read in infile

while (<>) {
	my @record = split;
	my $num_records = scalar(@record);
  # don't count as Hi-C contact if fails mapq or intra fragment test
	my $countme = 1;

	if (($record[1] eq $record[5]) && $record[3] == $record[7]) {
		$countme = 0;
	}
	elsif ($num_records > 8) {
		my $mapq_val = min($record[8],$record[11]);
		if ($mapq_val < $mapq_threshold) {
			$countme = 0;
		}
	}


	if ($countme) {
		my $frag_1 = &get_paired_fragment($record[0], $record[1], $record[2], $record[3], $record[9]);
		my $frag_2 = &get_paired_fragment($record[4], $record[5], $record[6], $record[7], $record[12]);
		print "${$frag_1}[0]\t${$frag_1}[1]\t${$frag_1}[2]\t${$frag_2}[0]\t${$frag_2}[1]\t${$frag_2}[2]\n";
	}
}







sub get_paired_fragment {
	if (!defined($chromosomes{$_[1]}[$index])) {next;}
	my $read_len; #set read length to adjust contact fragment ::: may need consider map cigar 
	my $cigar = $_[4];
	if($cigar =~ /(\d+)M/){
		$read_len = $1;
	}
	my $index = $_[3];
	my $five_prime_border;
	my $three_prime_border;
	if ($index == 0) {
		$five_prime_border = 1;
		$three_prime_border = $chromosomes{$_[1]}[$index];
	}else {
		$five_prime_border = $chromosomes{$_[1]}[$index-1];
		$three_prime_border = $chromosomes{$_[1]}[$index];
	}
	
	my @out;
	if ($_[0] == 0){
		@out = ($_[1], $_[2], $three_prime_border);
	}else{
		#avoid out boundary
		my $end = $_[2] + $read_len - 1 < $chromosomes{$_[1]}[-1] ?  $_[2] + $read_len - 1 : $chromosomes{$_[1]}[-1];
		@out = ($_[1], $five_prime_border, $end); #chromosome start end
	}
	my $out = \@out;
	return $out;
}
