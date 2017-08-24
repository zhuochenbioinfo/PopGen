use strict;
use warnings;

my($qFile,$indFile,$outFile) = @ARGV;

my $usage = "USAGE:\nperl $0 <Q file> <ind file> <output file>\n";
$usage .= "<Q file> is the output file of admixture with suffix Q.\n";
$usage .= "<ind file> is the individual list match with Q file by row.\n";

die $usage unless(@ARGV == 3);

# check the number of lines in qfile and indfile
my @lines1 = `wc -l $qFile|awk '{print \$1}'`;
my @lines2 = `wc -l $indFile|awk '{print \$1}'`;

chomp($lines1[0]);
chomp($lines2[0]);

unless($lines1[0] == $lines2[0]){
	die "<Q file> and <ind file> must have the same row number.\n";
}

open(IN1,"<$qFile") or die $!;
open(IN2,"<$indFile") or die $!;
open(OUT,">$outFile");

print OUT "#sample\tqType\tqScore\tallScores\n";

my $input2 = <IN2>;
while(<IN1>){
	chomp;
	my @scores = split/\s/;
	chomp($input2);
	my $rank;
	my $value = 0;
	for(my $i = 0; $i < @scores; $i++){
		if($scores[$i] > $value){
			$rank = $i + 1;
			$value = $scores[$i];
		}
	}
	print OUT "$input2\t$rank\t$value\t".join("|",@scores)."\n";
	$input2 = <IN2>;
}
close OUT;
close IN1;
close IN2;
