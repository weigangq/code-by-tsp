#!/usr/bin/env perl
#  Calculates pairwise codon distances (with variable ti/tv ratio)
#  between pairs of single-substitution codons
#  stop codons excluded
#  Purpose: test error-minimizing effects of a code in positional differences & ti/tv biases
use strict;
use Inline::Files; # necessary to embed multiple data sets
use Data::Dumper;
use Getopt::Std;
use List::Util qw(shuffle);
my %opts;
getopts('f:p:i:s:hb:q:w', \%opts); # titv
die "Usage: $0\n\t[-h]\t(help)\n\t[-f 'sgc|code_file|STDIN']\t(code)\n\t[-s 'pair (default)|fit|path|coord']\t(stats)\n\t[-p 'pol|hydro|vol|iso', default 'grantham']\t(aa prop)\n\t[-i (ti/tv, default 5)]\n\t[-b: begin codon ('TTT') -q: panelty for 1st (50); -w: panelty for 2nd (100)]\t(options for path)\n" if $opts{h} || (@ARGV == 0 && keys %opts == 0);
my (%dist_aa, %codon_table, %seen_aa);
my %transitions = ('A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C');
my $prop_suffix = $opts{p} || 'grantham';
my $titv = $opts{i} || 5; # tranition/transversion bias
my @raw_bases = qw(T C A G);
my $stat = $opts{s} || 'pair';
my $code_suffix = $opts{f} ? ($opts{f} eq 'sgc' ? 'sgc' : $opts{f}) : 'STDIN';
my $start = $opts{b} || 'TTT';
my $panelty_1st = $opts{q} || 50; # relative to 3rd position
my $panelty_2nd = $opts{w} || 100;
my %transitions = ('A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C');
my %transversions = ('A' => 'TC', 'G' => 'TC', 'C' => 'AG', 'T' => 'AG');
my %nt_pos = ('T' => 1, 'C' => 2,  'A' => 3,  'G' => 4);
my %pos_nt = (1 => 'T', 2 => 'C',  3 => 'A',  4 => 'G');

if ($code_suffix eq 'sgc') {
    while (<SGC>) {
        next unless /^[ATCG]/;
        chomp;
        my ($codon, $aa) = split;
        $codon_table{$codon} = $aa;
        $seen_aa{$aa}++;
    }
} elsif ($code_suffix eq 'STDIN') {
    while (<>) {
        next unless /^[ATCG]/;
        chomp;
        my ($codon, $aa) = split;
        $codon_table{$codon} = $aa;
        $seen_aa{$aa}++;
    }
} else {
   open IN, "<", $opts{f};
   while(<IN>) {
   next unless /^[ATCG]/;
   chomp;
   my ($codon, $aa) = split;
   $codon_table{$codon} = $aa;
   $seen_aa{$aa}++;
   }
   close IN;
}

while(<GRANTHAM>){
    chomp;
    next unless /^[A-Z]/;
    my ($aa1, $aa2, $dist) = split;
    $dist_aa{$aa1}{$aa2} = $dist;
}

my (%aa_vars, %aa_order);
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

my @aas = sort keys %seen_aa;
my @codons = sort keys %codon_table;
unless ($prop_suffix eq 'grantham') {
    for(my $i=0; $i<$#aas; $i++) {
	   for(my $j=$i+1; $j<=$#aas; $j++) {
	    $dist_aa{$aas[$i]}{$aas[$j]} = sprintf "%.2f", abs($aa_vars{$aas[$i]}->{$prop_suffix} - $aa_vars{$aas[$j]}->{$prop_suffix});
	    $dist_aa{$aas[$j]}{$aas[$i]} = sprintf "%.2f", abs($aa_vars{$aas[$i]}->{$prop_suffix} - $aa_vars{$aas[$j]}->{$prop_suffix});
	   }
    }    
    
    my %aa_val;
    foreach my $aa (@aas) {
        $aa_val{$aa} = $aa_vars{$aa}->{$prop_suffix};
    }    
    my @aas_sorted = sort {$aa_val{$a} <=> $aa_val{$b}} keys %aa_val;
    for (my $i=0; $i<=$#aas_sorted; $i++) {
        $aa_order{$aas_sorted[$i]} = $i+1;
    }
}

&print_code_fitness() if $stat eq 'fit'; # code fitness by pol, hydro, vol, or iso
&print_codon_dist() if $stat eq 'pair';
my @codon_path;
&print_codon_path() if $stat eq 'path';
&print_codon_coord() if $stat eq 'coord';

exit;

sub print_codon_coord { 
# to be filled
}

sub print_codon_path {
#    my @codon_path;
    &walk_codon_path(\$start);
    my $tour_length_codon = 0;
    my $tour_length_aa = 0;
    my @tour;
    my $last = $codon_path[-1];
    my $idx=1;
    foreach my $pair (@codon_path) {
        $tour_length_codon += &dist_codon($pair);
        $tour_length_aa += $dist_aa{$codon_table{$pair->{start}}}{$codon_table{$pair->{dest}}};
        push @tour, {   'codon'=>$pair->{start}, 
                        'aa' => $codon_table{$pair->{start}},
                        'aa_linear_idx' => $prop_suffix eq 'grantham' ? undef : $aa_order{$codon_table{$pair->{start}}},
                        'codon_d2'=>$pair->{type} eq 'i' ? 1:$titv, 
                        'codon_d1'=>$pair->{pos} == 2 ? $panelty_2nd : ($pair->{pos} == 1 ? $panelty_1st : 1),
                        'codon_linear_idx' => $idx++,
                    };
        push @tour, {   'codon'=>$pair->{dest}, 
                        'aa' => $codon_table{$pair->{dest}},
                        'aa_linear_idx' => $prop_suffix eq 'grantham' ? undef : $aa_order{$codon_table{$pair->{dest}}},
                        'codon_d2'=>0, 
                        'codon_d1'=>0,
                        'codon_linear_idx' => 0,
                    } if $pair eq $last;                    
#        $idx++;
#        push @tour, $pair->{start} . ";" . $pair->{type} . ";" . $pair->{pos};
#        push @tour, $pair->{dest} . ";" . "end" . ";" . "end" if $pair eq $last;       
    }
    foreach (@tour) { 
        print join ("\t", ($start, $_->{codon}, $_->{codon_d1}, $_->{codon_d2}, $_->{codon_linear_idx}, $_->{aa}), $_->{aa_linear_idx}), "\n";
    }
    warn $start, "\t", $tour_length_codon, "\t", $tour_length_aa, "\t", $prop_suffix, "\n";
#print Dumper(\@codon_path);
}

sub dist_codon {
    my $ref = shift;
    my $d = $ref->{type} eq 'i' ? 1 : $titv;
    $d += $panelty_2nd if $ref->{pos} == 2;
    $d += $panelty_1st if $ref->{pos} == 1;
    return $d;
}

sub walk_codon_path {
    my $ref = shift;
    my $cd = $$ref;
    my $last = &walk_1st($cd);
    push @codon_path, {'start' => $last, 'dest' => &sub_ti_2nd($last), 'type' => 'i', 'pos' => 2};
    my $last_from_2nd_ti1 = &walk_1st(&sub_ti_2nd($last)); # print "==>\n";
    push @codon_path, {'start' => $last_from_2nd_ti1, 'dest' => &sub_tv_2nd($last_from_2nd_ti1), 'type' => 'v', 'pos' => 2};
    my $last_from_2nd_tv = &walk_1st(&sub_tv_2nd($last_from_2nd_ti1)); # print "==>\n";
    push @codon_path, {'start' => $last_from_2nd_tv, 'dest' => &sub_tv_2nd($last_from_2nd_tv), 'type' => 'i', 'pos' => 2};
    my $last_from_2nd_ti2 = &walk_1st(&sub_ti_2nd($last_from_2nd_tv)); # print "==>\n";
}

sub walk_1st {
    my $cd = shift;
    my $last_from_3rd = &walk_3rd($cd); # print "-->\n";
    push @codon_path, {'start' => $last_from_3rd, 'dest' => &sub_ti_1st($last_from_3rd), 'type' => 'i', 'pos' => 1};
    my $last_from_1st_ti1 = &walk_3rd(&sub_ti_1st($last_from_3rd)); # print "==>\n";
    push @codon_path, {'start' => $last_from_1st_ti1, 'dest' => &sub_tv_1st($last_from_1st_ti1), 'type' => 'v', 'pos' => 1};
    my $last_from_1st_tv = &walk_3rd(&sub_tv_1st($last_from_1st_ti1)); # print "-->\n";
    push @codon_path, {'start' => $last_from_1st_tv, 'dest' => &sub_ti_1st($last_from_1st_tv), 'type' => 'i', 'pos' => 1};
    my $last_from_1st_ti2 = &walk_3rd(&sub_ti_1st($last_from_1st_tv)); # print "--->\n";
    return $last_from_1st_ti2;
}

sub sub_tv_2nd {
    my $cd = shift;
    my @bases = split //, $cd;
    my $second = $bases[1];
    my $second_tv = (split //, $transversions{$second})[0];
    return $bases[0] . $second_tv . $bases[2];
}

sub sub_ti_2nd {
    my $cd = shift;
    my @bases = split //, $cd;
    my $second = $bases[1];
    my $second_ti = $transitions{$second};
    return $bases[0] . $second_ti . $bases[2];
}

sub sub_tv_1st {
    my $cd = shift;
    my @bases = split //, $cd;
    my $last_two = $bases[1] . $bases[2];
    my $first = $bases[0];
    my $first_tv = (split //, $transversions{$first})[0];
    return $first_tv . $last_two;
}

sub sub_ti_1st {
    my $cd = shift;
    my @bases = split //, $cd;
    my $last_two = $bases[1] . $bases[2];
    my $first = $bases[0];
    my $first_ti = $transitions{$first};
    return $first_ti . $last_two;
}

sub walk_3rd {
    my $cd = shift;
#    print $cd, "->";
    my @bases = split //, $cd;
    my $first_two = $bases[0] . $bases[1];
    my $third = $bases[2];
    my $third_ti = $transitions{$third};
    my $tvs = $transversions{$third};
    my @tv_nts = split //, $tvs;
    my $cd_ti  = $first_two . $third_ti;
    my $cd_tv1 = $first_two . $tv_nts[0];
    my $cd_tv2 = $first_two . $tv_nts[1];
    push @codon_path, {'start' => $cd, 'dest' => $cd_ti, 'type' => 'i', 'pos' => 3};
    push @codon_path, {'start' => $cd_ti, 'dest' => $cd_tv1, 'type' => 'v', 'pos' => 3};
    push @codon_path, {'start' => $cd_tv1, 'dest' => $cd_tv2, 'type' => 'i', 'pos' => 3};
#    print $cd_ti, "=>", $cd_tv1, "->", $cd_tv2;
    return $cd_tv2;
}

sub print_code_fitness {
#    die "set {p} option to pol|hydro|iso|vol for codon fitness\n" if $prop_suffix eq 'grantham'; 
    my $n=0;
    my $codon_fitness = 0;
    foreach my $codon_i (@codons) {
	next if $codon_i eq 'X';
#        my $xi = $aa_vars{$codon_table{$codon_i}}->{$prop_suffix} unless $prop_suffix eq 'grantham';
        my $aai = $codon_table{$codon_i};
        my $nabes = &_get_neighbors($codon_i);
        foreach my $codon_j (@$nabes) {
	    next if $codon_j eq 'X';
#	       my $xj = $aa_vars{$codon_j->{'aa'}}->{$prop_suffix} unless $prop_suffix eq 'grantham';
	    my $aaj = $codon_j->{'aa'};
	    my $wt = $codon_j->{'ti'} ? 1 : $titv;
	    $codon_fitness += $wt * $dist_aa{$aai}{$aaj};
	    $n++;
        }
    }
    printf "%.4f\n", $codon_fitness/$n;
}

sub print_codon_dist {
    foreach (my $i = 0; $i < $#codons; $i++) {
        my @bases_x = split //, $codons[$i];
        my $aa_x = $codon_table{$codons[$i]};
        next if $aa_x eq 'X';
        foreach (my $j = $i+1; $j <= $#codons; $j++) {
	       my @bases_y = split //, $codons[$j];
	       my $aa_y = $codon_table{$codons[$j]};
	       next if $aa_y eq 'X';
	       my $sum_diff = 0;
	       my %diff;
	       for(my $k=1; $k<=3; $k++) {
	           next if $bases_x[$k-1] eq $bases_y[$k-1];
	           $diff{$k} = $transitions{$bases_x[$k-1]} eq $bases_y[$k-1] ? "ti" : "tv";
	           $sum_diff++;
	       }
           next if $sum_diff > 1;
	       my ($pos) = keys %diff;
	       my ($titv) = values %diff;
	       print join("\t", ($codons[$i], $codons[$j], $pos . "-". $titv, $dist_aa{$aa_x}{$aa_y}, $code_suffix, $prop_suffix)), "\n";
        }
    }
}

sub _get_neighbors {
    my $source = shift;
    my ($first, $second, $third) = split //, $source;
    my @neighbors;
    foreach my $rawbase (@raw_bases) {
        next if $rawbase eq $first;
        my $ti = ($transitions{$first} eq $rawbase) ? 1 : 0;
        push @neighbors, {
            'codon' => $rawbase . $second . $third,
            'ti' => $ti,
            'aa' => $codon_table{$rawbase . $second . $third}
        }
    }

    foreach my $rawbase (@raw_bases) {
        next if $rawbase eq $second;
        my $ti = ($transitions{$second} eq $rawbase) ? 1 : 0;
        push @neighbors, {
            'codon' => $first . $rawbase . $third,
            'ti' => $ti,
            'aa' => $codon_table{$first . $rawbase .  $third}
        }
    }

    foreach my $rawbase (@raw_bases) {
        next if $rawbase eq $third;
        my $ti = ($transitions{$third} eq $rawbase) ? 1 : 0;
        push @neighbors, {
            'codon' => $first . $second . $rawbase,
            'ti' => $ti,
            'aa' => $codon_table{$first . $second . $rawbase}
        }
    }
    return(\@neighbors);
}


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

# Standard genetic code ("X" for stop codon)
__SGC__
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
#########################################################
# Grantham Amino Acid chemical-property distance
#  modified by adding stop codon ("X"), with dist(X, non-X) = 500 (twice the max)
#  and dist(X,X) = 0
########################################################
__GRANTHAM__
A	A	0
A	C	195
A	D	126
A	E	107
A	F	113
A	G	60
A	H	86
A	I	94
A	K	106
A	L	96
A	M	84
A	N	111
A	P	27
A	Q	91
A	R	112
A	S	99
A	T	58
A	V	64
A	W	148
A	Y	112
C	A	195
C	C	0
C	D	154
C	E	170
C	F	205
C	G	159
C	H	174
C	I	198
C	K	202
C	L	198
C	M	196
C	N	139
C	P	169
C	Q	154
C	R	180
C	S	112
C	T	149
C	V	192
C	W	215
C	Y	194
D	A	126
D	C	154
D	D	0
D	E	45
D	F	177
D	G	94
D	H	81
D	I	168
D	K	101
D	L	172
D	M	160
D	N	23
D	P	108
D	Q	61
D	R	96
D	S	65
D	T	85
D	V	152
D	W	181
D	Y	160
E	A	107
E	C	170
E	D	45
E	E	0
E	F	140
E	G	98
E	H	40
E	I	134
E	K	56
E	L	138
E	M	126
E	N	42
E	P	93
E	Q	29
E	R	54
E	S	80
E	T	65
E	V	121
E	W	152
E	Y	122
F	A	113
F	C	205
F	D	177
F	E	140
F	F	0
F	G	153
F	H	100
F	I	21
F	K	102
F	L	22
F	M	28
F	N	158
F	P	114
F	Q	116
F	R	97
F	S	155
F	T	103
F	V	50
F	W	40
F	Y	22
G	A	60
G	C	159
G	D	94
G	E	98
G	F	153
G	G	0
G	H	98
G	I	135
G	K	127
G	L	138
G	M	127
G	N	80
G	P	42
G	Q	87
G	R	125
G	S	56
G	T	59
G	V	109
G	W	184
G	Y	147
H	A	86
H	C	174
H	D	81
H	E	40
H	F	100
H	G	98
H	H	0
H	I	94
H	K	32
H	L	99
H	M	87
H	N	68
H	P	77
H	Q	24
H	R	29
H	S	89
H	T	47
H	V	84
H	W	115
H	Y	83
I	A	94
I	C	198
I	D	168
I	E	134
I	F	21
I	G	135
I	H	94
I	I	0
I	K	102
I	L	5
I	M	10
I	N	149
I	P	95
I	Q	109
I	R	97
I	S	142
I	T	89
I	V	29
I	W	61
I	Y	33
K	A	106
K	C	202
K	D	101
K	E	56
K	F	102
K	G	127
K	H	32
K	I	102
K	K	0
K	L	107
K	M	95
K	N	94
K	P	103
K	Q	53
K	R	26
K	S	121
K	T	78
K	V	97
K	W	110
K	Y	85
L	A	96
L	C	198
L	D	172
L	E	138
L	F	22
L	G	138
L	H	99
L	I	5
L	K	107
L	L	0
L	M	15
L	N	153
L	P	98
L	Q	113
L	R	102
L	S	145
L	T	92
L	V	32
L	W	61
L	Y	36
M	A	84
M	C	196
M	D	160
M	E	126
M	F	28
M	G	127
M	H	87
M	I	10
M	K	95
M	L	15
M	M	0
M	N	142
M	P	87
M	Q	101
M	R	91
M	S	135
M	T	81
M	V	21
M	W	67
M	Y	36
N	A	111
N	C	139
N	D	23
N	E	42
N	F	158
N	G	80
N	H	68
N	I	149
N	K	94
N	L	153
N	M	142
N	N	0
N	P	91
N	Q	46
N	R	86
N	S	46
N	T	65
N	V	133
N	W	174
N	Y	143
P	A	27
P	C	169
P	D	108
P	E	93
P	F	114
P	G	42
P	H	77
P	I	95
P	K	103
P	L	98
P	M	87
P	N	91
P	P	0
P	Q	76
P	R	103
P	S	74
P	T	38
P	V	68
P	W	147
P	Y	110
Q	A	91
Q	C	154
Q	D	61
Q	E	29
Q	F	116
Q	G	87
Q	H	24
Q	I	109
Q	K	53
Q	L	113
Q	M	101
Q	N	46
Q	P	76
Q	Q	0
Q	R	43
Q	S	68
Q	T	42
Q	V	96
Q	W	130
Q	Y	99
R	A	112
R	C	180
R	D	96
R	E	54
R	F	97
R	G	125
R	H	29
R	I	97
R	K	26
R	L	102
R	M	91
R	N	86
R	P	103
R	Q	43
R	R	0
R	S	110
R	T	71
R	V	96
R	W	101
R	Y	77
S	A	99
S	C	112
S	D	65
S	E	80
S	F	155
S	G	56
S	H	89
S	I	142
S	K	121
S	L	145
S	M	135
S	N	46
S	P	74
S	Q	68
S	R	110
S	S	0
S	T	58
S	V	124
S	W	177
S	Y	144
T	A	58
T	C	149
T	D	85
T	E	65
T	F	103
T	G	59
T	H	47
T	I	89
T	K	78
T	L	92
T	M	81
T	N	65
T	P	38
T	Q	42
T	R	71
T	S	58
T	T	0
T	V	69
T	W	128
T	Y	92
V	A	64
V	C	192
V	D	152
V	E	121
V	F	50
V	G	109
V	H	84
V	I	29
V	K	97
V	L	32
V	M	21
V	N	133
V	P	68
V	Q	96
V	R	96
V	S	124
V	T	69
V	V	0
V	W	88
V	Y	55
W	A	148
W	C	215
W	D	181
W	E	152
W	F	40
W	G	184
W	H	115
W	I	61
W	K	110
W	L	61
W	M	67
W	N	174
W	P	147
W	Q	130
W	R	101
W	S	177
W	T	128
W	V	88
W	W	0
W	Y	37
Y	A	112
Y	C	194
Y	D	160
Y	E	122
Y	F	22
Y	G	147
Y	H	83
Y	I	33
Y	K	85
Y	L	36
Y	M	36
Y	N	143
Y	P	110
Y	Q	99
Y	R	77
Y	S	144
Y	T	92
Y	V	55
Y	W	37
Y	Y	0
A	X	500
X	A	500
C	X	500
X	C	500
D	X	500
X	D	500
E	X	500
X	E	500
F	X	500
X	F	500
G	X	500
X	G	500
H	X	500
X	H	500
I	X	500
X	I	500
K	X	500
X	K	500
L	X	500
X	L	500
M	X	500
X	M	500
N	X	500
X	N	500
P	X	500
X	P	500
Q	X	500
X	Q	500
R	X	500
X	R	500
S	X	500
X	S	500
T	X	500
X	T	500
V	X	500
X	V	500
W	X	500
X	W	500
Y	X	500
X	Y	500
X	X	0
