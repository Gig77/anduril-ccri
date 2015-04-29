#!/usr/bin/perl

# Positional Gene Enrichment (PGE) script
# Copyright (C) 2007 Roland Barriot
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#		
# email: 
#	roland.barriot@gmail.com
# mail: 
# 	Roland Barriot
# 	ESAT-SCD
# 	Katholieke Universiteit Leuven
# 	Kasteelpark Arenberg, 10home
# 	B-3001 Leuven
# 	Belgium

use Class::Struct;
use lib './';
use Hypergeometric;
#use Data::Dump qw(dump);


# PROGRAM PARAMETERS
struct ( Dataset => [
	name       => '$',
	loadMethod => '$',
	file       => '$',
]);

my $dataDir = '/users/sista/bioiuser/pge/data';
my %refDatasets;
# human
my %humanDatasets;
$refDatasets{'hgu95av2'}  = new Dataset(name=>'hgu95av2',loadMethod=>'loadAffy',file=>"$dataDir/HG_U95Av2.na22.annot.txt");
$humanDatasets{'hgu95av2'} = 1;
$refDatasets{'hgu133a'}   = new Dataset(name=>'hgu133a',loadMethod=>'loadAffy',file=>"$dataDir/HG-U133A.na22.annot.txt");
$humanDatasets{'hgu133a'} = 1;
$refDatasets{'hgu133a_subset'}   = new Dataset(name=>'hgu133a_subset',loadMethod=>'loadAffy',file=>"$dataDir/hgu133a_subset.txt");
$humanDatasets{'hgu133a_subset'} = 1;
$refDatasets{'ensembl42'} = new Dataset(name=>'ensembl42',loadMethod=>'load_ID_CHR_START_END_DESC',file=>"$dataDir/ensembl42.txt");
$humanDatasets{'ensembl42'} = 1;
$refDatasets{'entrez'} = new Dataset(name=>'entrez',loadMethod=>'load_ID_CHR_START_END_DESC',file=>"$dataDir/entrez.txt");
$humanDatasets{'entrez'} = 1;
$refDatasets{'symbols'} = new Dataset(name=>'symbols',loadMethod=>'load_ID_CHR_START_END_DESC',file=>"$dataDir/symbols.txt");
$humanDatasets{'symbols'} = 1;
$refDatasets{'refseq_dna'} = new Dataset(name=>'refseq_dna',loadMethod=>'load_ID_CHR_START_END_DESC',file=>"$dataDir/refseq_dna.txt");
$humanDatasets{'refseq_dna'} = 1;
$refDatasets{'refseq_peptide'} = new Dataset(name=>'refseq_peptide',loadMethod=>'load_ID_CHR_START_END_DESC', file=>"$dataDir/refseq_peptide.txt");
$humanDatasets{'refseq_peptide'} = 1;

# mouse
my %mouseDatasets;
$refDatasets{'mouse430_2'} = new Dataset(name=>'mouse430_2',loadMethod=>'loadAffy',file=>"$dataDir/Mouse430_2.na24.annot.csv");
$mouseDatasets{'mouse430_2'} = 1;

# user-defined
$refDatasets{'user'} = new Dataset(name=>'user',loadMethod=>'load_ID_CHR_START_END_DESC',file=>"__not_specified");
$humanDatasets{'user'} = 1;

my $dataset   	= 'ensembl42';
my $alpha 		= 0.05;
my $summarize 	= 3;
my $min_ratio   = 1.5;
my $mappedIDs 	= 'ensembl';

use vars qw/ %opt /;
use Getopt::Std;
my $opt_string = 'hlgr:q:f:a:s:c:m:';
getopts( "$opt_string", \%opt ) or usage();
usage() if $opt{h};

sub usage {
print STDERR << "EOF";

    usage: $0 [-hl] -q query_IDs [-r dataset] [-f dataset_file] [-a alpha_threshold] [-f redundancy_filter [-s ratio]] [-c chromosome] [-m mapped_IDs]

	-h			: this (help) message
	-l			: list available reference datasets
	-g			: allow gaps
	
	-q query_IDs_file
	
	-r reference dataset	: default $dataset
	-f str			: file name of user-defined reference dataset (required if -r 'user')
	-a float		: threshold for p-value significance. default: $alpha
	-s float		: pertinent if percentage*s > larger region percentage. default: $min_ratio
	-m str			: map to which IDs [ensembl|symbols]. symbols only available for Affy ref datasets. default: $mappedIDs
	-c int			: restrict search on particular chromosome. default: all
	 

    example: $0 -r ensembl42 -q my_file_with_ensembl_IDs

EOF
	exit;
}

# LIST AVAILABLE REFERENCE DATASETS
list_datasets() if ($opt{l}); 
sub list_datasets {
	print "Available datasets:\n\t- ".join("\n\t- ",keys %refDatasets);
	print "\n";
	exit;
}
#QUERY FILE
usage() if !$opt{q};
my $queryFile = $opt{q};

# REFERENCE DATASET
$dataset = lc($opt{r}) if $opt{r};
$dataset = lc("$dataset");
if (!defined $refDatasets{$dataset}) {
	print "Unknown reference dataset '$dataset'";
	list_datasets();
}
if ($dataset eq 'user') {
	if (!defined $opt{f}) {
		print STDERR "ERROR: User-defined reference data set not specified (-f)\n";
		exit;
	}
	if (! -e $opt{f}) {
		print STDERR "ERROR: User-defined reference data set ".$opt{f}." not found\n";
		exit;
	}
	$refDatasets{'user'} = new Dataset(name=>'user',loadMethod=>'load_ID_CHR_START_END_DESC',file=>$opt{f});
}

# ALPHA
$alpha = $opt{a} if $opt{a};
# REDUNDANCY FILTER
$min_ratio = $opt{s} if $summarize eq 3 && $opt{s};
# MAP TO WHICH IDs
$mappedIDs = $opt{m} if $opt{m};
# RESTRICT TO 1 CHROMOSOME
my $restrict_to_chromosome = $opt{c} if $opt{c};
# COMPUTATIONAL STUFF
my $di = Hypergeometric->new();
if ($opt{p} && $opt{p} eq 'mphyper') {
	$di = ModifiedHypergeometric->new();
}


# OUTPUT
my $err='';
my $unmapped = '';
my $raw='';
my $bed='';
my $query_symbols = '';


# DATA STRUCTURES
struct ( ID => [
	id     		=> '$',
	mappedID 	=> '$',
	symbol 		=> '$',
	ensg   		=> '$',
	chr    		=> '$',
	start 		=> '$',
	end    		=> '$',
	rank   		=> '$',
]);

struct ( Chromosome => [
	name  => '$',
	ids   => '@',
	size  => '$',
	query => '@',
	sEnriched => '@', # SIGNIFICANTLY ENRICHED
	dEnriched => '@',
	location  => '%', # GIVEN THE RANK/POSITION YOU GET A REFERENCE TO A GENES
]);
my %chrs; # CHROMOSOMES
if ($restrict_to_chromosome) {
	$chrs{$restrict_to_chromosome} = new Chromosome(name=>$restrict_to_chromosome, size=>0);
}
else {
	if (defined $humanDatasets{$dataset}) {
		for (my $i=1;$i<23;$i++) { $chrs{$i} = new Chromosome(name=>$i, size=>0); }
		$chrs{'x'} = new Chromosome(name=>'x', size=>0);
		$chrs{'y'} = new Chromosome(name=>'y', size=>0);
	}
	elsif (defined $mouseDatasets{$dataset}) {
		for (my $i=1;$i<20;$i++) { $chrs{$i} = new Chromosome(name=>$i, size=>0); }
		$chrs{'x'} = new Chromosome(name=>'x', size=>0);
		$chrs{'y'} = new Chromosome(name=>'y', size=>0);
	}
}

sub sortchr {
	if ($a->chr eq 'x' || $a->chr eq 'y' || $b->chr eq 'x' || $b->chr eq 'y') {
		return $a->chr cmp $b->chr;
	}
	return $a->chr <=> $b->chr;
}

struct ( Hit => [
	name      => '$',
	pvalue    => '$',
	pvalueadj => '$',
	chr       => '$',
	start     => '$',
	end       => '$',
	common    => '$',
	size      => '$',
	seqStart  => '$',
	seqEnd    => '$',
	pertinent => '$',
]);

# the real job
sub comparePertinentTargetSet {
	my $qset = shift; my @loc = @$qset;
	my $q = shift;
	my $popsize = shift;
	my $chr = shift;
	my @results;	
	my $qschr = scalar @loc; # SIZE OF THE QUERY ON THIS CHROMOSOME
	my $nchr = scalar @{$chrs{$chr}->ids};
	for (my $i=0;$i<$qschr-1;$i++) {
		# PERTINENCE TEST RULE 3
		next if ($i>0 &&  # NOT THE FIRST POSITION FOR LOWER BOUND
			$loc[$i-1]+1 eq $loc[$i]); # PREVIOUS I ITERATION WAS MORE PERTINENT
		for (my $j=$i+1;$j<$qschr;$j++) {
			# PERTINENCE TEST RULE 3
			next if ($j<$qschr-1 && # NOT THE LAST POSITION FOR UPPER BOUND
					$loc[$j]+1 eq $loc[$j+1]); # NEXT J ITERATION IS MORE PERTINENT
		
			my $t = $loc[$j]-$loc[$i]+1; # NUMBER OF ELEMENT IN THAT INTERVAL ON CHROMOSOME
			my $c = $j-$i+1; # NUMBER OF COMMON ELEMENTS
			my $pvalue = 1;
			if ($opt{p} && $opt{p} eq 'homemade') {
				$pvalue = 0;
				my $sum_inf = 0;
				my $sum_sup = 0;
				for (my $iq=$c;$iq<=$qschr;$iq++) {
					my $prod_inf = 1;
					my $prod_sup = 1;
					for (my $jq=0;$jq<=$c;$jq++) {
						$prod_inf *= ($nchr-$loc[$i]-$jq+1) / ($nchr-$jq);
						$prod_sup *= ($loc[$j]-$jq) / ($nchr-$jq);
					}
					$sum_inf += $prod_inf; if ($sum_inf>1) { $sum_inf=1; }
					$sum_sup += $prod_sup; if ($sum_sup>1) { $sum_sup=1; }
				}
				my $sum = 0;
				for (my $iq=$qschr;$iq<=$q && $iq <= $nchr;$iq++) {
					my $prod = 1;
					for (my $jq=0;$jq<$iq;$jq++) {
						$prod *= ($nchr-$jq)/($popsize-$jq);
					}
					$sum += $prod;
				}
				$pvalue = $sum_inf * $sum_sup * $sum;
			}
			else {
				$pvalue = $di->compute($popsize,$q,$t,$c);
			}
			my $r = new Hit(name=>$chr.','.$loc[$i].':'.$loc[$j], 
							pvalue=>$pvalue, 
							chr=>$chr, 
							start=>$loc[$i], 
							end=>$loc[$j],
							common=>$c,
							size=>$t,
							pertinent=>1,
							);
			push @results, $r;
		}
	}
	return @results;
}

# LOAD QUERY
my @queryElm;
open(F,"<$queryFile") || die ("Cannot open file for input '$queryFile'");
my $querySize=0;
while (<F>) {
	chomp;
	my @ar = split /\s/;
	foreach (@ar) {
		if (length $_>0) {
			push @queryElm, lc($_);		
			$querySize++;
		}
	}
}
close(F);

my %ids = &{$refDatasets{$dataset}->loadMethod}($refDatasets{$dataset}->file, $mappedIDs);

sub loadAffy {
	my $filename = shift;
	my $mappedIDs = shift;
	my %ids;
	my %alreadyMappedIDs;
	open(F,"<".$filename) || die ("cannot open ".$filename);
	while (<F>) {
		next if $. == 1;
		my @f = split '","';	
		my ($refID) = $f[0] =~ /.*"(\S+)/;
		$refID = lc($refID);
		my $loc = $f[12];
		next if $loc eq '---'; # NO LOCATION PROVIDED
		next if index($loc,"///") ne -1; # ONLY UNAMBIGUOUS PROBES
		my ($chr,$start,$end)     = $loc =~ /chr(\S+):(\d+)-(\d+)\s.*/;
		$chr = lc($chr);
		next if !defined $chrs{$chr}; # WE DON'T CONSIDER THIS CHROMOSOME
		
		my $ensg = $f[17];
		my $symbol = $f[14];
		if ($symbol eq '---') { $symbol = $ensg; }

		my $mappedID = $ensg;
		if ($mappedIDs eq 'symbols') {
			$mappedID = $symbol;
		}
		next if $mappedID eq '---' || index($mappedID,'///')>0; # KEEP ONLY PROBES THAT HAVE EXACTLY ONE MAPPED ID

		my $id = new ID(id=>$refID, mappedID=>$mappedID, 
						symbol=>$symbol, ensg=>$ensg, 
						chr=>$chr, start=>$start, end=>$end);
		push @{$chrs{$chr}->ids},$id unless defined $alreadyMappedIDs{$mappedID};
		$alreadyMappedIDs{$mappedID} = 1;
		
		my @trans;
		if (defined $ids{$refID}) {
			@trans = @{$ids{$refID}}; 
		}
		push @trans, $id;
		$ids{$refID} = \@trans;
	}
	close(F);
	return %ids;
}

sub load_ID_CHR_START_END_DESC {
	my $filename = shift;
	my %ids;
	open(F,"<".$filename) || die ("cannot open ".$filename);
	while (<F>) {		
		chomp;
		my ($refID, $chr, $start, $end, $desc) = split "\t";
		$refID = lc($refID);
		$chr = lc($chr);
		if (defined $chrs{$chr}) {
			my $id = new ID(id=>$refID, mappedID=>$refID, 
							symbol=>$refID, ensg=>$refID, 
							chr=>$chr, start=>$start, end=>$end);
			push @{$chrs{$chr}->ids},$id;
			my @trans;
			if (defined $ids{$refID}) {
				@trans = @{$ids{$refID}}; 
			}
			push @trans, $id;
			$ids{$refID} = \@trans;
		}
	}
	close(F);
	return %ids;
}


# COMPUTE RANKS OF GENES ON CHROMOSOMES
my %ids_ranks;
foreach my $c (values %chrs) {
	my $r=1;
	foreach my $g (sort { $a->start <=> $b->start } @{$c->ids}) {
		if (!defined $ids_ranks{$g->mappedID}) {
			$g->rank($r);
# 2015-04-27 | CF | accessing hash this way results in wrong region coordinates in script output! (not sure why)
#                   replaced with alternative syntax using accessor function
#			${%{$c->location}}{$r}=$g;
			$c->location($r, $g);
			$ids_ranks{$g->mappedID} = $r;
			$r++;
		}
		else {
			$g->rank($ids_ranks{$g->mappedID});
		}
	}
	$c->size($r-1);
}

# RETRIEVE QUERY ELEMENT LOCATIONS
my %queryElmHash;
foreach my $e (@queryElm) {
	if (!defined $ids{$e}) {
		$unmapped.="$e\n";
	}
	else {
		my $ai = $ids{$e};		
		foreach my $i (@$ai) {
			if (!defined $queryElmHash{$i->mappedID}) {
				push @{$chrs{$i->chr}->query}, $i;
				$queryElmHash{$i->mappedID} = 1;
				if (!defined $i->rank) {
					$i->rank($ids_ranks{$i->mappedID});
				}
			}
			$query_symbols.=$e.','.$i->mappedID.','.$i->chr.','.$i->start.','.$i->end.','.$i->ensg."\n";
		}
	}
}
$querySize = scalar keys %queryElmHash;

# PERFORM SEARCH
my @pvalues;
my @results; # UNFILTERED PERTINENT COMPARISON RESULTS (INCLUDING NON SIGNIFICANT SETS)
if ($opt{g}) { # ALLOW GAPPED REGIONS ( WITH NO GENES OF INTEREST OF SIZE LARGER THAN EXPECTED ON AVG)
	foreach my $c (values %chrs) {
		my @queryInt;
		my $lastRank=0;
		foreach	(sort {$a->rank <=> $b->rank} @{$c->query}) {
			if ($lastRank eq 0 || $lastRank ne $_->rank) {
				push @queryInt, $_->rank;
			}
			$lastRank = $_->rank;
		}	
		my @cresults = comparePertinentTargetSet(\@queryInt, $querySize, (scalar keys %ids), $c->name);
		foreach my $r (@cresults) {
			push @pvalues, $r->pvalue;
			push @results, $r;
		}
	}
}
else {
	my $avg_gap = (scalar keys %ids)/$querySize; # ONE GENE OF INTEREST EVERY ...
	foreach my $c (values %chrs) {
		my @queryInt = ();
		my @c_query = sort {$a->rank <=> $b->rank} @{$c->query};		
		next if scalar @c_query < 2;
		my $gene = shift @c_query;		
		push @queryInt, $gene->rank;
		my $lastRank=$gene->rank;		
		while (scalar @c_query>0) {
			my $gene = shift @c_query;
			if ($lastRank ne $gene->rank) {
				if ($gene->rank-$lastRank > $avg_gap) { # GAP HERE
					if (scalar @queryInt > 1) {
						my @cresults = comparePertinentTargetSet(\@queryInt, $querySize, (scalar keys %ids), $c->name);
						foreach my $r (@cresults) {
							push @pvalues, $r->pvalue;
							push @results, $r;
						}
					}
					@queryInt = ();				
				}
				push @queryInt, $gene->rank;
			}
			$lastRank = $gene->rank;
		}
		# LAST INTERVAL WITHOUT GAPS
		if (scalar @queryInt > 1) {
			my @cresults = comparePertinentTargetSet(\@queryInt, $querySize, (scalar keys %ids), $c->name);
			foreach my $r (@cresults) {
				push @pvalues, $r->pvalue;
				push @results, $r;
			}
		}
	}
}

# ADJUST SIGNIFICANCE
# FALSE DISCOVERY RATE
################
my $nbSignificant = 0;
sub min3 {
	my $m = shift;
	my $a = shift;
	my $b = shift;
	if ($a<$m) { $m = $a; }
	if ($b<$m) { $m = $b; }
	return $m;
}
sub FDR {
	my $alpha = shift;
	my $h     = shift;
	my @hits  = @$h;
	my @shits = sort {$a->pvalue <=> $b->pvalue} @hits;
	my $n = scalar @shits;
	if ($n > 0) {
		$shits[$n-1]->pvalueadj($shits[$n-1]->pvalue);
		for (my $i=$n-1;$i>=1;$i--) {
			$shits[$i-1]->pvalueadj( min3($shits[$i]->pvalueadj, $n/$i*$shits[$i-1]->pvalue, 1) );
		}	
	}
}

FDR($alpha,\@results);

# FILTER SIGNIFICANT HITS
my @sresults; # SIGNIFICANT PERTINENT HITS
foreach my $r (sort {$a->pvalueadj <=> $b->pvalueadj} @results) {
	if ($r->pvalueadj <= $alpha) {		
		push @sresults, $r;
		$nbSignificant++;
		push @{$chrs{$r->chr}->sEnriched}, $r;
	}
}

# FILTER REDUNDANCY (NESTED HITS WITH WORSE PVALUES)
foreach my $c (values %chrs) {
	my @hits = sort {$a->start <=> $b->start} @{$c->sEnriched}; # HITS SORTED BY START POSITION
	my $nb = scalar @hits;
	for (my $i=0;$i<$nb-1;$i++) {
		my $h1 = $hits[$i];
		for (my $j=$i+1;$j<$nb && $hits[$j]->start < $h1->end;$j++) {
			my $h2 = $hits[$j];				
			if ($h2->end <= $h1->end) { 								# H1 INCLUDES H2 I.E. H1 IS LARGER
				if ($h1->pvalue <= $h2->pvalue &&						# AND HAS A BETTER PVALUE
					$min_ratio*$h1->common/$h1->size > $h2->common/$h2->size) { # THEN KILL H2 IF RATIO TOO SMALL
					$h2->pertinent(0);
				}
				if ($h1->pvalue >= $h2->pvalue) { 	# AND H1 HAS AN EQUAL OR WORSE PVALUE => WORSE RATIO
					$h1->pertinent(0);
				}				
			}
			elsif ($h1->start == $h2->start) { # WE CANNOT USE THE SORT CRITERION, WE MUST TRY H1 INCLUDED IN H2
				if ($h1->end <= $h2->end) { 								# H1 INCLUDED IN H2 I.E. H1 IS LARGER
					if ($h2->pvalue <= $h1->pvalue &&						# AND HAS A BETTER PVALUE
						$min_ratio*$h2->common/$h2->size > $h1->common/$h1->size) { # THEN KILL H2 IF RATIO TOO SMALL
						$h1->pertinent(0);
					}
					if ($h2->pvalue >= $h1->pvalue) { 	# AND H1 HAS AN EQUAL OR WORSE PVALUE => WORSE RATIO
						$h2->pertinent(0);
					}				
				}				
			}
		}
	}
}


# PERFORMS OUTPUT
print "chr\tstart\tend\tpvalue\tpvalueadj\tcommon\tsize\tgenes\n";
$bed.= 	"track name=\"PGE\" ".
		"description=\"Positional Gene Enrichment\" ".
		"visibility=full color=200,100,0 altColor=0,100,200 priority=20\n";
foreach my $c (values %chrs) {
	my @aStart;
	my @aEnd;
	my $bedtmp = '';
	foreach my $r (sort {$a->pvalue <=> $b->pvalue} @{$c->sEnriched}) {
		next if (!$r->pertinent); # filter redundant hits
		
# 2015-04-27 | CF | accessing hash this way results in wrong region coordinates in script output! (not sure why)
#                   replaced with alternative syntax using accessor function
#		my $start = ${%{$c->location}}{$r->start}->start;
#		my $end   = ${%{$c->location}}{$r->end}->end;
		my $start = $c->location($r->start)->start;
		my $end   = $c->location($r->end)->end;
		my $score = 'inf';
		if ($r->pvalue > 0) {
			$score = int (-3 * log($r->pvalue)); # for the genome browser color
		}

		# find gene ids in enriched region
		my @common_gene_ids;
		for my $g (@{$c->query}) {
			if ($g->start >= $start && $g->end <= $end) {
				#print STDERR $g->id,"\n";
				push(@common_gene_ids, uc($g->id));
			}
		}
		
		print uc($c->name)."\t".$start."\t".$end."\t".$r->pvalue."\t".$r->pvalueadj."\t".$r->common."\t".$r->size."\t".join(",", @common_gene_ids)."\n";
		

#		$bedtmp.= 	"chr".uc($c->name)."\t".$start."\t".$end."\t".$r->pvalue."\n";
#		$raw.= uc($c->name)."\t".$start."\t".$end."\t".
#			$r->pvalue."\t".
#			$r->pvalueadj."\t".
#			$r->common."\t".
#			$r->size.
#			"\n";
	}
#	if (length $bedtmp > 0) {
#		$bed .= $bedtmp;
#	}
}

#print "<mapping>$query_symbols</mapping><raw>$raw</raw><bed>$bed</bed><alpha>$alpha</alpha><adjust>$adjustmentMethod</adjust><reference>$dataset</reference>";
