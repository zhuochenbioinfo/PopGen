use strict;
use warnings;
use Getopt::Long;

my($vcf,$outName,$cov,$maf,$pass,$pass1,$keepList,$regList);
my $usage = "\nTransform vcf to treemix SNP format.\nBy Zhuo Chen\nEmail: zhuochenbioinfo\@gmail.com\n\n";
$usage .= "Notice: the program picks only bi-allelic SNPs.\n";
$usage .= "USAGE:\nperl $0 --in <input vcf> --out <output file> --keep <keep sample list>\n";
$usage .= "<input vcf> is the input vcf file. [Necessary]\n";
$usage .= "<output file> is the name of output file. [Necessary]\n";
$usage .= "--keep <keep sample list>: each line consists of <sample name> and <group>. [Necessary]";
$usage .= "FILTERING OPTIONS:\n";
$usage .= "--cov <minium coverage ratio>, [0,1], Default=0\n";
$usage .= "--maf <minium allele frequency>, [0,1], Default=0\n";
$usage .= "--bed <keep region list> in BED format.\n";
$usage .= "--pass to keep SNPs with filter TAG [PASS] or [SnpCluster].\n";
$usage .= "--pass1 to keep SNPs with filter TAG [PASS].\n\n";

GetOptions(
	"in=s" => \$vcf,
	"out=s" => \$outName,
	"cov=s" => \$cov,
	"maf=s" => \$maf,
	"pass!" => \$pass,
	"pass1!" => \$pass1,
	"keep=s" => \$keepList,
	"bed=s" => \$regList,
) or die $usage;

die $usage unless(defined $vcf and defined $outName and defined $keepList);

unless(defined $cov){
	$cov = 0;
}
unless(defined $maf){
	$maf = 0;
}

my %hash_keep;
my %hash_group;
open(IN,"<$keepList") or die $!;
while(<IN>){
	chomp;
	my($sample,$group) = split/\t/;
	$hash_keep{$sample}{group} = $group;
	$hash_group{$group}{count} ++;
}
close IN;

my %hash_bed;
my %hash_chr; # This hash is used to check the end of the region list and stop reading vcf

if(defined $regList){
	open(IN,"<$regList") or die $!;
	while(<IN>){
		chomp;
		my($chr,$start,$end,$other) = split/\t/;
		next if(exists $hash_bed{$chr}{$start} and $hash_bed{$chr}{$start} > $end);
		$hash_bed{$chr}{$start} = $end;
	}
	close IN;
	# merge regions
	foreach my $chr(keys %hash_bed){
		$hash_chr{$chr} = "";
		my @starts = sort {$a <=> $b} keys %{$hash_bed{$chr}};
		for(my $i = 0; $i < @starts; $i++){
			my $starti = $starts[$i];
			my $endi = $hash_bed{$chr}{$starti};
			for(my $j = $i+1; $j < @starts; $j++){
				my $startj = $starts[$j];
				my $endj = $hash_bed{$chr}{$startj};
				last if($startj > $endi);
				if($endj > $endi){
					$endi = $endj;
					$hash_bed{$chr}{$starti} = $endj;
				}
				delete($hash_bed{$chr}{$startj});
				splice(@starts,$j,1);
				$j--;
			}
		}
	}
}


open(IN,"<$vcf") or die $!;
open(OUT,">$outName");

my @allSamples = ();
my @keepRanks = ();
my @remained_chrs = ();
my $chr_tmp = "";
my @regions = ();

while(<IN>){
	chomp;
	next if($_ =~ /^##/);
	if($_ =~ /^#CHROM/){
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@allSamples) = split/\t/;
		for(my $i = 0; $i < @allSamples; $i++){
			my $sample = $allSamples[$i];
			next unless(exists $hash_keep{$sample});
			my $group = $hash_keep{$sample}{group};
			push @{$hash_group{$group}{ranks}}, $i;
			push @keepRanks, $i;
		}
		foreach my $group(sort keys %hash_group){
			unless(exists $hash_group{$group}{ranks}){
				die "# ERROR: No sample in group:$group, please check your files!\n";
			}else{
				print "# CHECKING: Group:$group ".@{$hash_group{$group}{ranks}}." in $hash_group{$group}{count}\n";
			}
		}
		print OUT join(" ", sort keys %hash_group)."\n";
		next;
	}
	
	#last if($. == 5000);

	# split line
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$datas_join) = split/\t/,$_,10;

#=this function is under construction
	
	# Check if the position locates in the candidate regions
	goto NOBED unless(defined $regList);
	
	next unless(exists $hash_bed{$chr});
	if($chr ne $chr_tmp){
		@regions = ();
		foreach my $start(sort {$a <=> $b} keys %{$hash_bed{$chr}}){
			my $end = $hash_bed{$chr}{$start};	
			push @regions, "$start,$end";
		}
		$chr_tmp = $chr;
		if(exists $hash_chr{$chr}){
			delete($hash_chr{$chr});
		}
		@remained_chrs = keys %hash_chr;
		print "\t# Reading chr:$chr\n";
		if(@regions == 0){
			print "\t# Skipping chr:$chr\n";
		}
	}
	
	if(@remained_chrs == 0 and @regions == 0){
		last;
	}

	if(@regions == 0){
		next;
	}
	
	my $pickIt = 0;
	for(my $i = 0; $i < @regions; $i++){
		my($start,$end) = split/,/,$regions[$i];
		if($pos < $start){
			last;
		}elsif($pos > $end){
			splice(@regions, $i, 1);
			$i--;
			next;
		}else{
			$pickIt = 1;
			last;
		}
	}

	next unless($pickIt == 1);
	NOBED:

#=cut

	# Check if the variant is a bi-allelic SNP
	next unless(length($ref) == 1 and length($alt) == 1);
	
	# Check if the SNP
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	if(defined $pass1){
		next unless($filter eq "PASS");
	}
	
	# Check maf and coverage
	my @datas = split/\t/, $datas_join;

	my @keepDatas = ();
	foreach my $rank(@keepRanks){
		push @keepDatas, $datas[$rank];
	}
	my($cov_,$maf_) = get_cov_maf(@keepDatas);
	next unless($cov_ > $cov and $maf_ > $maf);
	
	# transform format
	my @values = ();
	foreach my $group(sort keys %hash_group){
		my $refcount = 0;
		my $altcount = 0;
		foreach my $rank(@{$hash_group{$group}{ranks}}){
			my $data = $datas[$rank];
			my $geno = vcf2geno($data);
			next if($geno == 9);
			$refcount += $geno;
			$altcount += (2 - $geno);
		}
		push @values, "$refcount,$altcount";
	}
	print OUT join(" ", @values)."\n";
}
close IN;
close OUT;

sub vcf2geno{
	my($data) = @_;
	my $geno = 9;
	if($data =~ /(\d+)\/(\d+)/){
		if($1 == $2 and $1 == 0){
			$geno = 2;
		}elsif($1 == $2 and $1 != 0){
			$geno = 0;
		}else{
			$geno = 1;
		}
	}
	return($geno);
}


sub get_cov_maf{
	my @genos = @_;
	my $altcount = 0;
	my $covcount = 0;
	my $lostcount = 0;
	
	foreach my $geno(@genos){
		my $gt = "-";
		if($geno =~ /^(\d+)\/(\d+)/){
			$gt = $1;
			if($1 ne $2){
				$gt = "h";
				$altcount ++;
			}elsif($1 != 0){
				$altcount += 2;
			}
		}
		if($gt ne '-'){
			$covcount += 2;
		}else{
			$lostcount += 2;
		}
	}
	# debug 20171117
	my $cov = 1 - ($lostcount/@genos)/2;
	my $maf = ($altcount/@genos)/2;
	my $maf2 = $cov - $maf;
	if($maf > $maf2){
		$maf = $maf2;
	}
	return($cov,$maf);
}
