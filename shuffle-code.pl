#!/usr/bin/env perl
# Different ways to shuffle/permute a genetic code (SGC or evolved)
use strict;
use Inline::Files; # necessary to embed multiple data sets
use Data::Dumper;
use Getopt::Std;
use List::Util qw(shuffle);
my %opts;
getopts('f:s:phve', \%opts); # file input or 'sgc', shuffle 1,2,3,4,5 (1-3: by position; 4=>21 blocks; 5=>random), assign aa (for evovled code) by polarity/hydropathy/volume/isoelectricity
die "Usage: $0 -f <sgc|evoled-code> (required) [ Optional: -s <1|2|3|4(blocks)|5(random)> -p(olarity; default, for assigning blocks to an evolved code (Traveling Codon Problem, TCP), not used for shuffling) -h(ydropathy) -v(olume) -e(iso-electricity)]\n" unless $opts{f};
my (%codon_table, %aa_blocks);
while(<GENCODE>) {
    next unless /^[ATCG]/;
    chomp;
    my ($codon, $aa) = split;
    $codon_table{$codon} = $aa;
    if ($aa_blocks{$aa}) { 
	my $ref = $aa_blocks{$aa};
	push @$ref, $codon;
    } else {
	$aa_blocks{$aa} = [ ($codon) ];
    }
}

my %aa_vars;
while(<AA>){
    chomp;
    next unless /^[A-Z]/;
    my @data = split;
    $aa_vars{$data[0]} = {
	'pol' => $data[1],
	'hydro' => $data[2],
	'vol' => $data[3],
	'iso' => $data[4]
    };
}

# the following are necessary for assigning aa to evolved code
my $aa_p = 'pol'; # default polarity
$aa_p = 'pol' if $opts{p};
$aa_p = 'hydro' if $opts{h};
$aa_p = 'vol' if $opts{v};
$aa_p = 'iso' if $opts{e};

my @aa_sorted = sort {$aa_vars{$a}->{$aa_p} <=> $aa_vars{$b}->{$aa_p}} keys %aa_vars;

unless ($opts{f} eq 'sgc') { # input an evolved code
    open INCODE, "<", $opts{f};
    %aa_blocks = ();
    my %cd_pos;
    while(<INCODE>) {
	next unless /[ATCG][ATCG][ATCG]/;
	chomp;
	my @data = split;
	my ($codon, $pos) = ($data[1], $data[0]);
	$cd_pos{$codon} = $pos;
    }

    foreach my $codon (sort keys %cd_pos) {
	my $aa = $aa_sorted[$cd_pos{$codon}-1]; # ordered by aa traits
	$codon_table{$codon} = $aa;
	if ($aa_blocks{$aa}) { 
	    my $ref = $aa_blocks{$aa};
	    push @$ref, $codon;
	} else {
	    $aa_blocks{$aa} = [ ($codon) ];
	}
    }
    close INCODE;
}

&shuffle_all() if $opts{s} == 5;
&shuffle_aa_blocks() if $opts{s} == 4;
&shuffle_by_pos($opts{s}) if $opts{s} <= 3 && $opts{s} >= 1;

foreach my $cd (sort keys %codon_table) {
    print $cd, "\t", $codon_table{$cd}, "\n";
}
#print Dumper(\%codon_table); exit;

exit;

sub shuffle_by_pos {
    my $pos = shift;
    my %codon_table_degenerated;
    foreach my $cd (keys %codon_table) {
	my @bases = split //, $cd;
	my $aa = $codon_table{$cd};
#	next if $aa eq 'X';
	my $cd_deg;
	for (my $i=0; $i<=2; $i++) {
	    $cd_deg .= $bases[$i] unless $pos == $i+1;
	}
	if ($codon_table_degenerated{$cd_deg}) {
	    my $ref = $codon_table_degenerated{$cd_deg};
	    push @$ref, $aa;
	} else { 
	    $codon_table_degenerated{$cd_deg} = [ ($aa) ];
	}
    }
#    print Dumper(\%codon_table_degenerated);

    my @bicodes = shuffle(keys %codon_table_degenerated);
    my @aas = shuffle(values %codon_table_degenerated);
    my %shuffled_blks;
    for (my $i=0; $i<=$#bicodes; $i++) {
	$shuffled_blks{$bicodes[$i]} = $aas[$i];
    }

    foreach my $bi (sort keys %shuffled_blks) {
	my @cd_3Letters;
	my @bases = split //, $bi;
	foreach my $new_base (qw(A T C G)) { # restore 3-letter codes from 2-letter codes
	    if ($pos == 1) {
		push @cd_3Letters, $new_base . $bi;
	    } elsif ($pos == 2) {
		push @cd_3Letters, $bases[0] . $new_base . $bases[1];
	    } else {
		push @cd_3Letters, $bi . $new_base;
	    }
	}

	my @aa = shuffle (@{$shuffled_blks{$bi}}); # randomize

	for (my $i=0; $i<=3; $i++) {
	    $codon_table{$cd_3Letters[$i]} = $aa[$i];
	}
    }
#    return \%codon_table;
}

sub shuffle_all {
    my $ref = shift;
#    my %table = %$ref;
    my @aas = shuffle(values %codon_table);
    my @codons = sort keys %codon_table;
    for (my $i=0; $i<=$#codons; $i++) {
	$codon_table{$codons[$i]} = $aas[$i];
    }
#    return $ref;
}

sub shuffle_aa_blocks {
    my @aa = shuffle(keys %aa_blocks);
    my @blk = shuffle(values %aa_blocks);
    my %shuffled_blks;
    for (my $i=0; $i<=$#aa; $i++) {
	$shuffled_blks{$aa[$i]} = $blk[$i];
    }

    foreach my $aa (sort keys %shuffled_blks) {
	foreach my $cd (@{$shuffled_blks{$aa}}) {
	    $codon_table{$cd} = $aa;
	}
    }
#    return \%codon_table;
#    print Dumper(\%shuffled_blks);
}

# Standard genetic code ("X" for stop codon)
__GENCODE__
GGG	G
GGA	G
GGT	G
GGC	G
GAG	E
GAA	E
GAT	D
GAC	D
GTG	V
GTA	V
GTT	V
GTC	V
GCG	A
GCA	A
GCT	A
GCC	A
AGG	R
AGA	R
AGT	S
AGC	S
AAG	K
AAA	K
AAT	N
AAC	N
ATG	M
ATA	I
ATT	I
ATC	I
ACG	T
ACA	T
ACT	T
ACC	T
TGG	W
TGA	X
TGT	C
TGC	C
TAG	X
TAA	X
TAT	Y
TAC	Y
TTG	L
TTA	L
TTT	F
TTC	F
TCG	S
TCA	S
TCT	S
TCC	S
CGG	R
CGA	R
CGT	R
CGC	R
CAG	Q
CAA	Q
CAT	H
CAC	H
CTG	L
CTA	L
CTT	L
CTC	L
CCG	P
CCA	P
CCT	P
CCC	P

##########################################################
#  Woese's AA polarity & other AA one-dimensional measures
#  From Haig & Hurst (1991). Table 1
#  https://link.springer.com/article/10.1007/BF02103132
# Polarity, Hydropathy, volume, isoelectricity
########################################################
__AA__
A	7.0 1.8 31 6.00
R	9.1 -4.5 124 10.76
D	13.0 -3.5 54 2.77
N	10.0 -3.5 56 5.41
C	4.8 2.5 55 5.07
E	12.5 -3.5 8.3 3.22
Q	8.6 -3.5 85 5.65
G	7.9 -0.4 3 5.97
H	8.4 -3.2 96 7.59
I	4.9 4.5 111 6.02
L	4.9 3.8 111 5.98
K	10.1 -3.0 119 9.74
M	5.3 1.9 105 5.74
F	5.0 2.8 132 5.48
P	6.6 -1.6 32.5 6.30
S	7.5 -0.8 32 5.68
T	6.6 -0.7 61 6.16
W	5.2 -0.9 170 5.89
Y	5.4 -1.3 136 5.66
V	5.6 4.2 84 5.96
X	0 0 0
