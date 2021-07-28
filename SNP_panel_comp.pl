#!perl;
use strict; 
#use warnings; 
use v5.8.8;
use Data::Dumper;

## Script for the comparison between LGC genotyping file and VCF for Pengelly et al., 2013 SNP panel.
## Completed 2Dec14, Reuben Pengelly. 

## Open input files.
if ($ARGV[0] !~ m/.vcf$/) {die "ERROR: $ARGV[0] does not appear to be a VCF file.\n";}
open VCF, $ARGV[0] or die "ERROR: Unable to open VCF file \'$ARGV[0]\'.\n";
open GENO, $ARGV[1] or die "ERROR: Unable to open LGC genotype file \'$ARGV[1]\'.\n";


## Load SNP details into @rray and %ash for easy rsID comparison.
## A %ash from each perspective (location and rsID) is made for ready access.
chomp(my @panel = <DATA>);
my %panelLoc; my %panelRs; my @panelLoc;
foreach (@panel) {
	my @splitLine = split /\t/, $_;
	my $locus = join ':', @splitLine[0,1];
	$panelLoc{$locus} = $splitLine[2];
	$panelRs{$splitLine[2]} = $locus;
	push @panelLoc, $locus;
}


## Split VCF into header (discarded), field names and data sections.
my (@vcfSamples, @vcfGenoLines, %vcfGeno);
while (<VCF>) {
	chomp $_;
	if ($_ =~ m/^##/) {next;}
	if (($_ =~ m/^#/) && ($_ !~ m/^##/)) {
		my @head = split /\t/, $_;
		@vcfSamples = @head[9..$#head];
	}
	if ($_ !~ m/#/) {
		push @vcfGenoLines, $_;
	}
}
close VCF;


## Parse VCF genotypes into HoH.
my $j=0; ## $j used to ensure HoH initalised once only.
foreach my $line (@vcfGenoLines) {
	my @splitLine = split /\t/, $line;
	my $i = 0;
	my $loc = join ':', @splitLine[0,1];
	if (!exists $panelLoc{$loc}) { ## Skip VCF genotypes not in panel.
		print STDERR "\nWARNING: \'$loc\' in VCF file is not defined in panel, skipping.\n";
		next;
	}
	foreach (@splitLine[9..$#splitLine]) { ## Replace VCF coding with CSV allele list.
		my $id = $vcfSamples[$i]; 
		$_ =~ m/([01\.][|\/][01\.])/;  		
		my $geno = $1;
		$geno =~ s/0/$splitLine[3]/g; ## 1 & 0 coding replaced by alleles as defined in VCF line.
		$geno =~ s/1/$splitLine[4]/g;
		$geno =~ s/[\.23]/N/g; ## Ignore missing and non bi-allelic alleles, unlikely to occur aside by errors.
		my @genos = split /\/|/, $geno;
		@genos = sort @genos;
		$geno = join ',', @genos;
		
		if ($j == 0) { ## Parse into HoH
			$vcfGeno{$id} = { $loc => $geno }; ## Initialise anonymous HoH;
			$j++;
		} else {
			$vcfGeno{$id}{$loc} = $geno; ## Add to anonymous HoH.
		}
		$i++;
	}
}


## Count VCF samples.
my $vcfSampleN = @vcfSamples;
my $vcfGenoN = @vcfGenoLines;
print STDOUT "Running on $vcfGenoN loci for $vcfSampleN samples from VCF file.\n";


## Parse LGC genotypes and split.
my @genoRs; my %genoGeno; my $headered = 0; my $samplei = 0; my $assayed = 0;
while (<GENO>) {
	chomp $_;
	$_ =~ s/\r//g;
	if (($. == 1) && ($_ !~ 'KBiosciences grid report')) {
		die "\nERROR: Genotype file does not appear to be of KBiosciences format, check input.\n\n";
	}
	if ($headered == 0) { ## Parse header with genotype names if not seen yet. Skip all other header lines.
		if ($_ =~ /DNA \\ Assay/) {
			$headered = 1;
			@genoRs = split /,/, $_;
			shift @genoRs;
			foreach (@genoRs) {
				if ($_ =~ m/^$/) {next;}
				if (!exists $panelRs{$_}) {
					print STDERR "WARNING: \'$_\' in genotypes does not apear in the panel, skipping.\n";
					next;
				} else {
					$assayed++;
					next;
				}
			}
		}
		next;	
	}
	$_ =~ s/\?/N:N/g;
	my @splitLine = split /,/, $_;
	my $id = shift @splitLine;

	## Sort alleles for consistency.
	map {my @splitGeno = split /:/, $_; 
		@splitGeno = sort @splitGeno; 
		$_ = join ',', @splitGeno;
	} @splitLine; 	

	## Parse into HoH.
	my $genoi = 0;
	foreach (@splitLine) {
		my $loc = $genoRs[$genoi];
		if (!exists $panelRs{$loc}) {
			next;
		} else {
			$loc = $panelRs{$loc};
		}
		if ($samplei == 0) { 			
			$genoGeno{$id} = { $loc => $splitLine[$genoi] }; ## Initialise anonymous HoH;
		} else {
			$genoGeno{$id}{$loc} = $splitLine[$genoi]; ## Add to anonymous HoH.
		}
		$genoi++;
	}
	$samplei++;
}
close GENO;


## Count genotyping samples.
my $genoGenoN = @genoRs;
$assayed = $assayed * 2;
print STDOUT "Running on $genoGenoN loci for $samplei samples from genotype file.\n";


## Check which samples are paired between data sources.
## Unpaired data are ignored downstream.
my $vcfSkipSample = 0;
foreach my $id (keys %vcfGeno) {
	if (!exists $genoGeno{$id}) {
		print STDERR "WARNING: No matched-sample found for VCF sample \'$id\', skipping.\n";
		delete $vcfGeno{$id};
		$vcfSkipSample++;
	}
}
print "Skipping $vcfSkipSample VCF samples as they have no pair.\n";
my $genoSkipSample = 0;
foreach my $id (keys %genoGeno) {
	if (!exists $vcfGeno{$id}) {
		print STDERR "WARNING: No matched-sample found for genotyping sample \'$id\', skipping.\n";
		delete $genoGeno{$id};
		$genoSkipSample++;
	}
}
print "Skipping $genoSkipSample genotyping samples as they have no pair.\n";


## Compare Genotypes and see if match; if not, flip non matching; if still not, print error and die once all SNPs checked .
my $error = 0; my $iter = 0;
foreach my $loc (keys %panelLoc) {
	my %alleles;
	my $attempt = 0;
	foreach my $id (keys %vcfGeno) { ## Make %ash with alleles as keys. 
		my $genoV = $vcfGeno{$id}{$loc};
		$alleles{$genoV} = $id . 'vcf';
		my $genoG = $genoGeno{$id}{$loc};
		$alleles{$genoG} = $id . 'geno';
	}
	#print $loc; print Dumper \%alleles; For debugging if unexpected alleles not matching error.
	$iter++;
	delete $alleles{'N,N'};
	my $alleleCount = keys %alleles;
	print Dumper %alleles;
	if ($alleleCount == 3) {
		next;
	} elsif (($alleleCount == 6) && ($attempt == 0)) {
		map {$genoGeno{$_}{$loc} =~ tr/ACGT/TGCA/; 
		my @splitGeno = split /,/, $genoGeno{$_}{$loc}; 
		@splitGeno = sort @splitGeno; 
		$genoGeno{$_}{$loc} = join ',', @splitGeno;} keys %genoGeno;
		$attempt++;
		redo;
	} else {
		print STDERR "ERROR: Alleles do not match, and are not complementary, between genotype and VCF for \'$loc\'.\n";
		$error++;
	}
}
if ($error != 0) {die;}


## Compare profile along individuals
open REPORT, ">SNP_comparisons.txt";
print REPORT "#SampleID	Alleles_assayed Alleles_available	Alleles_mismatching	Proportion_matching	Paired?\n";
open FASTA, ">SNP_comparisons.fasta";
$, = ','; print FASTA "## Panel SNP order: @panelLoc\n"; $, = '';

my @ids = keys %vcfGeno;
@ids = sort @ids;
foreach my $id (@ids) {
	my @vcfPrint; my @genoPrint;
	foreach my $locus (@panelLoc) {
		push @vcfPrint, $vcfGeno{$id}{$locus};
		push @genoPrint, $genoGeno{$id}{$locus};
	}
	# print "$id" . 'VCF  ' . "@vcfPrint\n";
	# print "$id" . 'geno ' . "@genoPrint\n";	

	my $vcfToPrint = join '', @vcfPrint;
	my $genoToPrint = join '', @genoPrint;
	$vcfToPrint =~ s/,//g; 
	$genoToPrint =~ s/,//g;
	print FASTA ">$id" . "_vcf\n" . "$vcfToPrint\n";
	print FASTA ">$id" . "_geno\n" . "$genoToPrint\n";

	my $missing = 0; my $mismatch = 0;
	for (my $i = 0; $i <= $#vcfPrint; $i++) {
		if (($vcfPrint[$i] =~ m/N/) || ($genoPrint[$i] =~ m/N/)) {
			$missing += 2;
			next;
		}
		my @vcfDip = split /,/, $vcfPrint[$i];
		if ($genoPrint[$i] !~ "$vcfDip[0],") {
			$mismatch++;
		}
		if ($genoPrint[$i] !~ ",$vcfDip[1]") {
			$mismatch++;
		}
	}
	
	my $prop;
	if ($assayed - $missing == 0) { ## Prevent div by 0 error if no available genotypes.
		$prop = 0;
	} else {
		$prop = (1 - ($mismatch / ($assayed - $missing)));
	}
	my $match = $assayed - $mismatch;
	my $state;
	if ($missing > 16) {
		$state = 'Insufficient_data';
	} elsif  ($prop == 1) {
		$state = 'Exact_match';
	} elsif (($prop > 0.9) && ($mismatch <= 2)) {
		$state = 'Approx_match';
	} else {
		$state = 'Non_matching';
	}
	my $available = $assayed - $missing;
	print REPORT "$id	$assayed	$available	$mismatch	$prop	$state\n";
}

close FASTA; close REPORT;


__DATA__ ## Taken from paper; do not edit.
1	179520506	rs1410592	NPHS2
1	67861520	rs2229546	IL12RB2
2	169789016	rs497692	ABCB11
2	227896976	rs10203363	COL4A4
3	4403767	rs2819561	SUMF1
4	5749904	rs4688963	EVC
5	82834630	rs309557	VCAN
6	146755140	rs2942	GRM1
7	48450157	rs17548783	ABCA13
8	94935937	rs4735258	PDP1
9	100190780	rs1381532	TDRD7
10	100219314	rs10883099	HPSE2
11	16133413	rs4617548	SOX6
12	993930	rs7300444	WNK1
13	39433606	rs9532292	FREM2
14	50769717	rs2297995	L2HGDH
15	34528948	rs4577050	SLC12A6
16	70303580	rs2070203	AARS
17	71197748	rs1037256	COG1
18	21413869	rs9962023	LAMA3
19	10267077	rs2228611	DNMT1
20	6100088	rs10373	FERMT1
21	44323590	rs4148973	NDUFV3
22	21141300	rs4675	SERPIND1