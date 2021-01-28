=credit
By Zhuo CHEN, IGDB, CAS
Contact: chenomics@163.com
zhuochen@genetics.ac.cn
=cut


use strict;
use warnings;
use Getopt::Long;

my($vcf,$outName,$cov,$maf,$pass,$pass1,$keepList,$regList,$renameSample,$stdchr,$chrNameList);
my $usage = "\nTransform vcf to tped format.\nBy Zhuo Chen\nEmail: chenomics\@163.com\n\n";
$usage .= "Notice: the program picks only bi-allelic SNPs.\n";
$usage .= "USAGE:\nperl $0 --in <input vcf> --out <output prefix>\n";
$usage .= "<input vcf> is the input vcf file. [Necessary]\n";
$usage .= "<output prefix> is the prefix of output files. [Necessary]\n";
$usage .= "FILTERING OPTIONS:\n";
$usage .= "--cov <minium coverage ratio>, [0,1], Default=0\n";
$usage .= "--maf <minium allele frequency>, [0,1], Default=0\n";
$usage .= "--keep <keep sample list>\n";
$usage .= "--bed <keep region list> in BED format.\n";
$usage .= "--pass to keep SNPs with filter TAG [PASS] or [SnpCluster].\n";
$usage .= "--pass1 to keep SNPs with filter TAG [PASS].\n";
$usage .= "\n--rs: set the option to rename sample for R.\n\n";
$usage .= "--stdchr: set the option if you need to re-code the chrName by keeping the numbers in the chrName. For example: Chr5 -> 5\n\n";

GetOptions(
	"in=s" => \$vcf,
	"out=s" => \$outName,
	"cov=s" => \$cov,
	"maf=s" => \$maf,
	"pass!" => \$pass,
	"pass1!" => \$pass1,
	"keep=s" => \$keepList,
	"bed=s" => \$regList,
	"rs!" => \$renameSample,
	"stdchr!" => \$stdchr,
	"chrlist=s" => \$chrNameList,
) or die $usage;

die $usage unless(defined $vcf and defined $outName);

unless(defined $cov){
	$cov = 0;
}
unless(defined $maf){
	$maf = 0;
}

my %hash_cnl;
if(defined $chrNameList){
	open(IN,"<$chrNameList") or die $!;
	while(<IN>){
		chomp;
		my($chr,$name) = split/\t/;
		$hash_cnl{$chr}{name} = $name;
	}
	close IN;
}

my %hash_keep;

if(defined $keepList){
	open(IN,"<$keepList") or die $!;
	while(<IN>){
		chomp;
		my($sample,$other) = split/\t/;
		$hash_keep{$sample} = "";
	}
	close IN;
}

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


if($vcf =~ /gz$|gzip$/){
	open(IN,"zcat $vcf|") or die $!;
}else{
	open(IN,"<$vcf") or die $!;
}
open(TFAM,">$outName.tfam");
open(TPED,">$outName.tped");
open(MAP,">$outName.map");
open(SNP,">$outName.snp");
open(CHR,">$outName.chr");

my @allSamples = ();
my @keepSamples = ();
my @keepRanks = ();
my @remained_chrs = ();
my $chr_tmp = "";
my @regions = ();
my $chrNum = 0;
my %hash_tchr;

while(<IN>){
	chomp;
	next if($_ =~ /^##/);
	if($_ =~ /^#CHROM/){
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@allSamples) = split/\t/;
		if(defined $keepList){
			for(my $i = 0; $i < @allSamples; $i++){
				my $sample = $allSamples[$i];
				next unless(exists $hash_keep{$sample});
				push @keepRanks, $i;
				if(defined $renameSample){
					$sample = renameSampleR($sample);
				}
				push @keepSamples, $sample;
			}
		}else{
			@keepRanks = (0..@allSamples-1);
			for(my $i = 0; $i < @allSamples; $i++){
				my $sample = $allSamples[$i];
				if(defined $renameSample){
					$sample = renameSampleR($sample);
				}
				push @keepSamples, $sample;
			}
		}
		if(@keepRanks == 0){
			die "#ERROR: No sample remain after filtering, please check your files!\n";
		}else{
			print "# Remain ".@keepRanks." samples after filtering.\n";
		}
		foreach my $sample(@keepSamples){
			print TFAM "$sample $sample 0 0 0 -9\n";
		}
		close TFAM;
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
	if(defined $pass1){
		next unless($filter eq "PASS");
	}
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	
	# Check maf and coverage
	my @datas = split/\t/, $datas_join;

	my @keepDatas = ();
	foreach my $rank(@keepRanks){
		push @keepDatas, $datas[$rank];
	}
	my($cov_,$maf_,$major) = get_cov_maf(@keepDatas);
	next unless($cov_ > $cov and $maf_ > $maf);
	
	# rename chr
	unless(exists $hash_tchr{$chr}){
		if(exists $hash_cnl{$chr}){
			$chrNum = $hash_cnl{$chr}{name};
		}elsif(defined $stdchr){
			($chrNum) = $chr =~ /(\d+)/;
		}else{
			$chrNum++;
		}
		$hash_tchr{$chr}{num} = $chrNum;
		print CHR "$chr\t$chrNum\n";
	}
	
	# transform format
	my $genoline = "$chrNum $chrNum:$pos 0 $pos";
	print MAP "$genoline\n";
	print SNP "$chr\t$pos\t$ref\t$alt\t$cov_\t$maf_\t$major\n";
	foreach my $rank(@keepRanks){
		my $data = $datas[$rank];
		my $geno = vcf2tped($data,$major);
		$genoline .= " ".$geno;
	}
	print TPED "$genoline\n";	
}
close IN;
close TPED;
close MAP;
close SNP;

sub vcf2geno{
	my($data) = @_;
	my $geno = "NA";
	if($data =~ /(\d+)\/(\d+)/){
		if($1 == $2 and $1 == 0){
			$geno = 1;
		}elsif($1 == $2 and $1 != 0){
			$geno = -1;
		}else{
			$geno = 0;
		}
	}
	return($geno);
}

sub vcf2tped{
	my($data,$major) = @_;
	die "#ERROR: major allele shall be coded as 0 or 1.\n" unless($major == 0 or $major == 1);
	my $geno = "0 0"; # missing data
	if($data =~ /(\d+)\/(\d+)/){
		if($1 != $2){
			$geno = "1 2";
		}else{
			$geno = "1 1";
			if($1 == $major){
				$geno = "2 2";
			}
		}
	}elsif($data =~ /(\d+)\|(\d+)/){
		if($1 != $2){
			$geno = "1 2";
		}else{
			$geno = "1 1";
			if($1 == $major){
				$geno = "2 2";
			}
		}
	}
	return($geno);
}

sub get_cov_maf{
	my @genos = @_;
	my $altcount = 0;
	my $covcount = 0;
	my $lostcount = 0;
	my $major = 0;
	
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
		}elsif($geno =~ /^(\d+)\|(\d+)/){
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
		$major = 1;
	}
	return($cov,$maf,$major);
}

sub renameSampleR{
	my($sample) = @_;
	if($sample =~ /^\d/){
		$sample = "V$sample";
	}
	$sample =~ s/\-/\_/g;
	return($sample);
}
