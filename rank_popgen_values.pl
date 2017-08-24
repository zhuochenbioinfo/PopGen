use strict;
use warnings;
use Getopt::Long;

my($in,$colrank,$out);
my $usage = "USAGE:\nperl $0 --in <input file> --rank <1-based column rank> --out <output file>\n";
$usage .= "<input file> is a tab-delimited file of windowed pogene values. The first row shall be colnames. [Necessary]\n";
$usage .= "<1-based column rank> is the colrank of the value used to rank the windows. Default is the last column. [Optional]\n";

GetOptions(
	"in=s" => \$in,
	"out=s" => \$out,
	"rank=s" => \$colrank,
) or die $usage;

die $usage unless(defined $in and defined $out);
unless(defined $colrank){
	$colrank = 0;
}

my @lines = ();
my %hash_value;
my $header;
$colrank--;

open(IN,"<$in") or die $!;
while(<IN>){
	chomp;
	if($. == 1){
		$header = $_;
		next;
	}
	my $rank = $. - 2;
	push @lines, $_;
	my @items = split/\t/;
	my $value = $items[$colrank];
	my @check = checkNumeric($value);
	next unless(@check == 1);
	$hash_value{$rank}{value} = $value;
}
close IN;

my @sortedRanks = sort {$hash_value{$a}{value} <=> $hash_value{$b}{value}} keys %hash_value;

for(my $i = 0; $i < @sortedRanks; $i++){
	my $score = ($i+1)/@sortedRanks;
	my $rank = $sortedRanks[$i];
	$score = sprintf("%.4f", $score);
	$hash_value{$rank}{score} = $score;
}

open(OUT,">$out");
print OUT "$header\tpercentage\n";
for(my $i = 0; $i < @lines; $i++){
	my $line = $lines[$i];
	my $score = "NA";
	if(exists $hash_value{$i}{score}){
		$score = $hash_value{$i}{score};
	}
	print OUT "$line\t$score\n";
}
close OUT;

sub checkNumeric{
	my @nums = @_;
	my @out = ();
	my $reg1 = qr/^-?\d+(\.\d+)?$/;
	my $reg2 = qr/^-?0(\d+)?$/;
	foreach my $num(@nums){
		next unless($num =~ $reg1 && $num !~ $reg2);
		push @out,$num;
	}
	return(@out);
}
