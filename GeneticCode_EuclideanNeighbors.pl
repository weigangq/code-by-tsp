#!/usr/bin/env perl
# Oliver: this code has embedded data sets for Gramtham & Standard Genetic code
#  So you don't need to read any external files
#  Also, there are functions to calculate codon and aa distances
use lib("/Users/oliverattie/Downloads/Switch-2.17/blib/lib/");
use lib("/Users/oliverattie/Downloads/Inline-Files-0.69/blib/lib/");

use Switch;
use Inline::Files; # necessary to embed multiple data sets
use Data::Dumper;
use Getopt::Std;
#my $exe, $trial, $iteration, $flag;
#my $start_time, $elapse_time, $new_E, $t, $T0, $T;
my %Boltz;
my $sum_B;
my %d;
my %distance;
my %V;
my %city;
my @visit;
my @tsp_len;
my @tsp_time;
#my $sum_len, $sum_time, $min_E, $max_E, $temp;

my %opts;
getopts('i:w:', \%opts);
my (%dist_codon, %dist_aa, %dist_combo, %codon_table);
my $titv = $opts{i} || 4; # transition:transversion ratio; default 4
my $wt_codon = $opts{w} || 10; # tune this up to increase the weight for codon distances
my %transitions = ('A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C');
while(<GENCODE>) {
    next unless /^[ATCG]/;
    chomp;
    my ($codon, $aa) = split;
    $codon_table{$codon} = $aa;
    $amino_acids1{$aa}=1;
}
#my @amino_acids=sort keys %amino_acids1;
#print Dumper(\%codon_table); exit;
open(CODON_FILE, "/Users/oliverattie/codon-path.txt");
my @codons1;
while(<CODON_FILE>){
    chomp;
    push(@codons1,$_);
}
close(CODON_FILE);
for(my $i=0; $i<@codons1; $i++){
    $CodonNum{$codons1[$i]}=$i;
}
while(<GRANTHAM>){
    chomp;
    next unless /^[A-Z]/;
    my ($aa1, $aa2, $dist) = split;
    $dist_aa{$aa1}{$aa2} = $dist;
}
my @amino_acids;
open(AA_SCALED, "/Users/oliverattie/aa-scaled.txt");
$header=<AA_SCALED>;
while(<AA_SCALED>){
    my @aa_data=split;
    $Polarity{$aa_data[0]}=$aa_data[1];
    $Hydrophobicity{$aa_data[0]}=0*$aa_data[2];
    $Volume{$aa_data[0]}=0*$aa_data[3];
    $Iso_electricity{$aa_data[0]}=$aa_data[4];
    push(@amino_acids, $aa_data[0]);
}
close(AA_SCALED);
my %d_aa;
for(my $aa1=0; $aa1<@amino_acids; $aa1++){
    for(my $aa2=0; $aa2<@amino_acids; $aa2++){
	$d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]}=sqrt(($Polarity{$amino_acids[$aa1]}-$Polarity{$amino_acids[$aa2]})*($Polarity{$amino_acids[$aa1]}-$Polarity{$amino_acids[$aa2]})+($Hydrophobicity{$amino_acids[$aa1]}-$Hydrophobicity{$amino_acids[$aa2]})*($Hydrophobocity{$amino_acids[$aa1]}-$Hydrophobocity{$amino_acids[$aa2]})+($Volume{$amino_acids[$aa1]}-$Volume{$amino_acids[$aa2]})*($Volume{$amino_acids[$aa1]}-$Volume{$amino_acids[$aa2]}));
    }
}
my $maximum=0;
for(my $aa1=0; $aa1<@amino_acids; $aa1++){
    for(my $aa2=0; $aa2<@amino_acids; $aa2++){
	if($d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]}>$maximum){
	    $maximum=$d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]};
	}
    }
}
for(my $aa1=0; $aa1<@amino_acids; $aa1++){
    for(my $aa2=0; $aa2<@amino_acids; $aa2++){
	$d_aa1{$aa1}{$aa2}=$d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]}/$maximum;
    }
}

#print Dumper(\%dist_aa);
# Codon distances: could be simplied by 1/2, but okay since only done once
#my @codons = sort keys %codon_table;
foreach my $codon_x (@codons1) {
    foreach my $codon_y (@codons1) {
	$dist_codon{$codon_x}{$codon_y} = &get_codon_dist($codon_x, $codon_y);
	$dist_combo{$codon_x}{$codon_y} = &get_combo_dist($codon_x, $codon_y, \%codon_table);
    }
}
my $trial=20;
my $iteration=200;
my $B=7;
my $D=2;
my $N1=64;
my $N=21;
#srand(10);

my $maximum=0;
#my @codons=sort keys %codon_table;
for(my $codons1=0; $codons1<$N1; $codons1++){
    for(my $codons2=0; $codons2<$N1; $codons2++){
#	print $dist_combo{$codons[$codon1]}{$codons[$codon2]}." ".$codons[$codon1]." ".$codons[$codon2]." ".$codon_table{$codons[$codon1]}." ".$codon_table{$codons[$codon2]}."\n";
	if($dist_codon{$codons1[$codons1]}{$codons1[$codons2]}>$maximum){
	$maximum=$dist_codon{$codons1[$codons1]}{$codons1[$codons2]};
	
	}
    }
}
my $maximum1=0;
for(my $codon1=0; $codon1<$N; $codon1++){
    for(my $codon2=0; $codon2<$N; $codon2++){
	if($dist_combo{$codons1[$codon1]}{$codons1[$codon2]}>$maximum1){
	    $maximum1=$dist_combo{$codons1[$codon1]}{$codons1[$codon2]};
	}
    }
}
#my $d_P=2*$maximum+10;
#my $B=$d_P;
for(my $aa1=0; $aa1<$N; $aa1++){
    for(my $aa2=0; $aa2<$N; $aa2++){
	$d{$aa1}{$aa2}=$dist_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]}/$maximum;
    }
}
for(my $codon1=0; $codon1<$N1; $codon1++){
    for(my $codon2=0; $codon2<$N1; $codon2++){
#	if($codon_table{$codons1[$codon1]} eq $codon_table{$codons1[$codon2]}){
#	    $d_codon{$codon1}{$codon2}=-3/4;
#	}else{
	$d_codon{$codon1}{$codon2}=$dist_combo{$codons1[$codon1]}{$codons1[$codon2]}/$maximum1;
#	}
    }
}
#    $d_Q=0;
#    for(my $codon2=0; $codon2<$N; $codon2++){
#	$d_Q+=$d{0}{$codon2}/($N-1);
#    }
#	$T_X=$d_P*$d_Q/$d{0}{1};
#	print $T_X." ";
#    print "\n";


#for($exe=0; $exe<$trial; $exe++){
#    srand(1);
    for(my $i=0; $i<$N; $i++){
	for(my $j=0; $j<$N1; $j++){
	    $V{$i}{$j}=0.1*rand(1);
	}
    }
my $T0=0;
#for(my $i=1; $i<=75; $i++){
$T=$T0=0.6;
#my $T=$T0=0.6;
#print "T:".$i."\n";
#    my $T=$T0=sqrt($d_P*$d_Q/$N);
#		$T=$T0=1.6;
    while($T>=0.4){
		my $t=0;
#		$T=$T0=0.6;

		$T=$T0=$T0/(2**0.33);
		for(my $k=0; $k<$iteration; $k++){
		    $t=$t+0.005;
		    for(my $x=0; $x<$N; $x++){
			$sum_B=0;
			for(my $i=0; $i<$N1; $i++){
			$Boltz{$x}{$i}=exp(-del_energy(\%d_aa1,\%V,$x,$i)/$T);
			$sum_B=$sum_B+$Boltz{$x}{$i};
			}
		    
		    for(my $i=0; $i<$N1; $i++){
			if($sum_B!=0){
			$V{$x}{$i}=$Boltz{$x}{$i}/$sum_B;
			}
#			print "V:".$V{$x}{$i}." ";
			
		    }
#			print "\n";
		}
		
		    my @visit;
		    my %output;
		    my $len=tour(\%d_aa1,\%V, \@visit,\%output);
#		    print "Tour Length:".$len."\n";

#		    print "Stop City X-Coord Y-Coord\n";
#		    for(my $i=0; $i<$N1; $i++){
#			my @data=split(/ /, $output{$i});
#			my $codon_data="";
#			for(my $j=0; $j<@data; $j++){
#			    $codon_data.=$codons1[$data[$j]]." ";
#			}
#			$neighborref=get_neighbors($codons1[$i]);
#			$n=@{$neighborref};
#			print $n."\n";
#			print $i." ".$codons1[$i]." ".$codons1[$visit[$i]]."  ".$visit[$i]." ".$codon_table{$codons1[$visit[$i]]}." ".$codon_data."\n";
#		    }
		$T=$T0/(1.0+$t);
		    print STDERR $T."\n";
		}

}
for(my $i=0; $i<$N; $i++){
    for(my $j=0; $j<$N1; $j++){
	if($V{$i}{$j}>=0.1){
	print "V:".$i." ".$j." ".$V{$i}{$j}." ";
	}
    }
    print "\n";
}
my @visit;
my %output;
my $len=tour(\%d_aa1, \%V, \@visit,\%output);
print "Tour Length:".$len*$maximum."\n";
print "Stop Codon Nucleotides Amino Acids\n";
for(my $i=0; $i<$N1; $i++){
#    my @data=split(/ /, $output{$i});
#    my $codon_data="";
#    for(my $j=0; $j<@data; $j++){
#    $codon_data.=$codons1[$data[$j]]." ";
#    }
    print $i." ".$codons1[$i]." ".$visit[$i]." ".$amino_acids[$visit[$i]]." ".$codon_data."\n";
}
#}
#sub transv_dist{
#    $nuc1=shift;
#    $nuc2=shift;
#    $dist=0;
#    if($nuc1 eq $nuc2){
#	$dist=0;
#    }
#    if((($nuc1 eq "A")&&($nuc2 eq "G"))||(($nuc1 eq "G")&&($nuc2 eq "A"))){
#	$dist=1;
#    }
#    if((($nuc1 eq "C")&&($nuc2 eq "T"))||(($nuc1 eq "T")&&($nuc2 eq "C"))){
#	$dist=1;
#    }
#    if((($nuc1 eq "A")&&($nuc2 eq "C"))||(($nuc1 eq "C")&&($nuc2 eq "A"))){
#	$dist=4;
#    }
#    if((($nuc1 eq "G")&&($nuc2 eq "T"))||(($nuc1 eq "T")&&($nuc2 eq "G"))){
#	$dist=4;
#    }
#    if((($nuc1 eq "A")&&($nuc2 eq "T"))||(($nuc1 eq "T")&&($nuc2 eq "A"))){
#	$dist=4;
#    }
#    if((($nuc1 eq "C")&&($nuc2 eq "G"))||(($nuc1 eq "G")&&($nuc2 eq "C"))){
#	$dist=4;
#    }
#    return $dist;
#}
#sub energy{
#    $dref=shift;
#    $Vref=shift;
#    my $term1=0;
#    my $term2=0;
#    for(my $i=0; $i<$N; $i++){
#	for(my $x=0; $x<$N; $x++){
#	    for(my $y=0; $y<$N; $y++){
#		if($y!=$x){
#		    $term1=$term1+$$Vref{$x}{$i}*$$Vref{$y}{$i};
#		}
#	    }
#	}
#    }
#	for(my $x=0; $x<$N; $x++){
#	    for(my $y=0; $y<$N; $y++){
#		if($y!=$x){
#		    for(my $i=0; $i<$N; $i++){
#			switch($i){
#			    case 0{$term2=$term2+$$dref{$x}{$y}*$$Vref{$x}{$i}*($$Vref{$y}{$i+1}+$$Vref{$y}{$N-1});break;}
#			    case ($N-1){$term2=$term2+$$dref{$x}{$y}*$$Vref{$x}{$i}*($$Vref{$y}{0}+$$Vref{$y}{$i-1});break;}
#			  else {$term2=$term2+$$dref{$x}{$y}*$$Vref{$X}{$i}*($$Vref{$y}{$i+1}+$$Vref{$y}{$i-1});break;}

#			}
#}
#}
#}
#}
#			return ($B*0.5*$term1+$D*0.5*$term2);
#}		    
	
sub del_energy{
    my $dref=shift;
    my $Vref=shift;
    $x=shift;
    $i=shift;
    my $term1=0;
    my $term2=0;
    for(my $y=0; $y<$N; $y++){
	if($y!=$x){
	    $term1=$term1+$$Vref{$y}{$i};
	}
    }
    for(my $y=0; $y<$N; $y++){
	if($y!=$x){
	    my $neighborref=get_neighbors($codons1[$i]);
	    for(my $k=0; $k<@{$neighborref}; $k++){
		$term2=$term2+($$dref{$x}{$y})*($$Vref{$y}{$CodonNum{$$neighborref[$k]}});

	    }
	}
    }
   return($B*$term1+$D*$term2);
}
sub tour{
    my $dref=shift;
    my $Vref=shift;
    my $visitref=shift;
    my $outputref=shift;
    my $len=0;
    my $ord=0;
    for(my $x=0; $x<$N; $x++){
	for(my $i=0; $i<$N1; $i++){
	    if($$Vref{$x}{$i}>=0.1){
		$$visitref[$i]=$x;
		$$outputref{$i}.=$x." ";
		$ord++;
	    }
	}
    }
    if($ord!=$N){
	return 0;
    }
    for(my $i=0; $i<$N-1; $i++){
	$len=$len+$$dref{$$visitref[$i]}{$$visitref[$i+1]};
    }
    $len=$len+$$dref{$$visitref[$N-1]}{$$visitref[0]};
    return $len;
}
sub get_neighbors{
    my $source=shift;
    my ($first,$second,$third)=split //, $source;
    my @neighbors;
    my @raw_bases=qw(A T C G);
    my %transition=("A"=>"G", "G"=>"A", "C"=>"T", "T"=>"C");
#    foreach my $rawbase(@raw_bases){
#	next if $rawbase eq $first;
#	my $ti=($transition{$first} eq $rawbase)?1:0;
#	if($ti){
#	push (@neighbors,$rawbase.$second.$third);
#	    'codon'=>$rawbase.$second.$third,
#	    'ti'=>$ti,
#	}
#    }
#    foreach my $rawbase(@raw_bases){
#        next if $rawbase eq $second;
#	my $ti=($transition{$second} eq $rawbase)?1:0;
#	if($ti){
#        push (@neighbors,$first.$rawbase.$third);
#	    'codon'=>$first.$rawbase.$third,
#	    'ti'=>$ti,
#	}
#    }
    foreach my $rawbase(@raw_bases){
	next if $rawbase eq $third;
	my $ti=($transition{$third} eq $rawbase)?1:0;
#	print $rawbase." ".$ti."\n";
#	if($ti){
	    push(@neighbors, $first.$second.$rawbase);
#	}
    }

    return \@neighbors;
}
sub get_combo_dist {
    my ($x, $y, $ref) = @_;
    my $dist_codon_pairs = $dist_codon{$x}{$y}; # this is too small (overwhelmed) compared with Gramtham distance below; perhaps add a weight to equlize the two effects (mutation disttance & AA distance)?
    my $dist_aa_pairs = $dist_aa{$ref->{$x}}{$ref->{$y}};
#    return [ ($dist_codon_pairs, $dist_aa_pairs) ];
    return $wt_codon * $dist_codon_pairs + $dist_aa_pairs;
}

sub get_codon_dist {
    my ($x, $y) = @_;
    return 0 if $x eq $y;
    my @bases_x = split //, $x;
    my @bases_y = split //, $y;
    my $sum = 0;
    for(my $i=0; $i<=2; $i++) {
	next if $bases_x[$i] eq $bases_y[$i];
	if($i==0){
	$sum += 25*$transitions{$bases_x[$i]} eq $bases_y[$i] ? 1 : $titv;
	}if($i==1){
	    $sum+=50*$transitions{$bases_x[$i]} eq $bases_y[$i] ? 1 : $titv;
	}if($i==2){
	    $sum+=$transitions{$bases_x[$i]} eq $bases_y[$i] ? 1 : $titv;
	}
    }
    return $sum;
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
