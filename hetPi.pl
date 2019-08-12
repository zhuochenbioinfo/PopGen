use strict;
use warnings;
use Getopt::Long;

my($vcf,$minmaf,$mincov,$out,$keepList,$window,$pass,$pass1);

my $usage = "USAGE:\nperl $0 --in <vcf> --out <output file>\n";
$usage .= "--keep <keep sample list>: [Optional]\n";
$usage .= "--maf: minor allele frequency [0,1], default=0 [Optional]\n";
$usage .= "--cov: minium coverage [0,1], default=0 [Optional]\n";
$usage .= "--window <window size>: default=10000 [Optional]\n";
$usage .= "--pass to keep SNPs with filter TAG [PASS] or [SnpCluster].\n";
$usage .= "--pass1 to keep SNPs with filter TAG [PASS].\n";

GetOptions(
	"in=s" => \$vcf,
	"out=s" => \$out,
	"cov=s" => \$mincov,
	"maf=s" => \$minmaf,
	"window=s" => \$window,
	"keep=s" => \$keepList,
	"pass!" => \$pass,
	"pass1!" => \$pass1,
) or die $usage;

die $usage unless(defined $vcf and defined $out);

unless(defined $mincov){
	$mincov = 0;
}
unless(defined $minmaf){
	$minmaf = 0;
}
unless(defined $window){
	$window = 10 * 1000;
}

my %hash_keep;
if(defined $keepList){
	open(IN,"<$keepList") or die $!;
	while(<IN>){
		chomp;
		my($sample,$other) = split/\t/;
		$hash_keep{$sample}{rank} = "";
	}
	close IN;
}

open(IN,"<$vcf") or die $!;
open(OUT,">$out");
print OUT "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\n";

my @pickedRanks;
my %hash_window;

my $tmp = "";

while(<IN>){
	chomp;
	next if($_ =~ /^##/);
	
	my($chr,$pos,$id,$ref,$alts_join,$qual,$filter,$info,$format,$datas_join) = split/\t/,$_,10;
	
	# dissect samples
	if($chr =~ /#CHROM/){
		my @samples = split/\t/,$datas_join;
		my $tail = @samples - 1;
		my @pickedSamples = ();
		
		print @pickedRanks."\n";
		unless(defined $keepList){
			@pickedRanks = (0..$tail);
			next;
		}
		for(my $i = 0; $i < @samples; $i++){
			next unless(exists $hash_keep{$samples[$i]});
			push @pickedRanks, $i;
		}
		NOKEEP:
		print "#Keeping ".@pickedRanks." samples from ".@samples." samples.\n";
		next;
	}
	
	my $wRank = int($pos/$window);
	my $wStart = $wRank * $window + 1;
	my $wEnd = $wRank * $window + $window;
	
	if($tmp ne "$chr|$wRank"){
		goto EMPTY if($tmp eq "");
		my($chr_tmp,$wRank_tmp) = $tmp =~ /(\S+)\|(\d+)/;
		goto EMPTY unless(exists $hash_window{$chr_tmp}{$wRank_tmp});
		my $wStart_tmp = $wRank_tmp * $window + 1;
		my $wEnd_tmp = $wRank_tmp * $window + $window;
		my $windowPi = $hash_window{$chr_tmp}{$wRank_tmp}{sum}/$window;
		my $count = @{$hash_window{$chr_tmp}{$wRank_tmp}{pi}};
		print OUT "$chr_tmp\t$wStart_tmp\t$wEnd_tmp\t$count\t$windowPi\n";
		EMPTY:
		$tmp = "$chr|$wRank";
	}
	
	
	# Check if the variant is a bi-allelic SNP
	next unless(length($ref) == 1 and length($alts_join) == 1);
	
	# Check if the SNP pass the filter
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	if(defined $pass1){
		next unless($filter eq "PASS");
	}
	
	my @datas = split/\t/,$datas_join;
	
	my @pickedDatas = ();
	foreach my $rank(@pickedRanks){
		push @pickedDatas, $datas[$rank];
	}
	my($cov,$maf,$pi) = hetCovMafPi(@pickedDatas);
	next unless($cov >= $mincov and $maf >= $minmaf);
	
	push @{$hash_window{$chr}{$wRank}{pi}}, $pi;
	$hash_window{$chr}{$wRank}{sum} += $pi;
}
close IN;

# output the last bin if exists
my($chr_tmp,$wRank_tmp) = $tmp =~ /(\S+)\|(\d+)/;
unless(exists $hash_window{$chr_tmp}{$wRank_tmp}){
	my $wStart_tmp = $wRank_tmp * $window + 1;
	my $wEnd_tmp = $wRank_tmp * $window + $window;
	my $windowPi = $hash_window{$chr_tmp}{$wRank_tmp}{sum}/$window;
	my $count = @{$hash_window{$chr_tmp}{$wRank_tmp}{pi}};
	print OUT "$chr_tmp\t$wStart_tmp\t$wEnd_tmp\t$count\t$windowPi\n";
}
close OUT;

sub hetCovMafPi{
	my @genos = @_;
	my $refcount = 0;
	my $altcount = 0;
	my $hetcount = 0;
	my $covcount = 0;
	
	foreach my $geno(@genos){
		if($geno =~ /^(\d+)\/(\d+)/){
			$covcount++;
			if($1 != $2){
				$hetcount++;
			}elsif($1 == 0){
				$refcount++;
			}else{
				$altcount++;
			}
		}
	}
	
	my $cov = $covcount/@genos;
	my @counts = sort {$b <=> $a} ($refcount,$hetcount,$altcount);
	my $maf = ($covcount - $counts[0])/@genos;
	my $refratio = 0;
	my $hetratio = 0;
	my $altratio = 0;
	if($covcount > 0){
		$refratio = $refcount/$covcount;
		$hetratio = $hetcount/$covcount;
		$altratio = $altcount/$covcount;
	}
	my $pi = $refratio * $hetratio * 1 + $altratio * $hetratio * 1 + $refratio * $altratio * 2;
	return($cov,$maf,$pi);
}
