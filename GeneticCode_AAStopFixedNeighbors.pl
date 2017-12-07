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
my %Boltz;
my $sum_B;
my %d;
my %distance;
my %V;
my %city;
my @visit;
my @tsp_len;
my @tsp_time;


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
}

my @codons1;
while(<CODON_FILE>){
    chomp;
    if($_ ne "TAG" && $_ ne "TAA" && $_ ne "TGA"){
    push(@codons1,$_);
    }
}
for(my $i=0; $i<@codons1; $i++){
    $CodonNum{$codons1[$i]}=$i;
}
my @amino_acids;
$header=<AA_SCALED>;
while(<AA_SCALED>){
    my @aa_data=split;
    $Polarity{$aa_data[0]}=$aa_data[1];
    $Hydropathy{$aa_data[0]}=$aa_data[2];
    $Volume{$aa_data[0]}=$aa_data[3];
    $Iso_electricity{$aa_data[0]}=$aa_data[4];
    push(@amino_acids, $aa_data[0]);

}
$Polarity{"X"}=-5;
$Hydropathy{"X"}=5;
$Volume{"X"}=-5;
my %d_aa;
for(my $aa1=0; $aa1<@amino_acids; $aa1++){
    for(my $aa2=0; $aa2<@amino_acids; $aa2++){
	$d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]}=sqrt(($Polarity{$amino_acids[$aa1]}-$Polarity{$amino_acids[$aa2]})*($Polarity{$amino_acids[$aa1]}-$Polarity{$amino_acids[$aa2]})+($Hydropathy{$amino_acids[$aa1]}-$Hydropathy{$amino_acids[$aa2]})*($Hydropathy{$amino_acids[$aa1]}-$Hydropathy{$amino_acids[$aa2]})+($Volume{$amino_acids[$aa1]}-$Volume{$amino_acids[$aa2]})*($Volume{$amino_acids[$aa1]}-$Volume{$amino_acids[$aa2]}));
    }
}
while(<AA_AS_CITIES>){
    chomp;
    push(@aa_as_cities, $_);
}
my $maximum=0;
for(my $aa1=0; $aa1<@amino_acids; $aa1++){
    for(my $aa2=0; $aa2<@amino_acids; $aa2++){
	if($d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]}>$maximum){
	    $maximum=$d_aa{$amino_acids[$aa1]}{$amino_acids[$aa2]};
	}
    }
}
for(my $aa1=0; $aa1<@aa_as_cities; $aa1++){
    for(my $aa2=0; $aa2<@aa_as_cities; $aa2++){
	$d{$aa1}{$aa2}=$d_aa{$codon_table{$codons1[$aa1]}}{$codon_table{$codons1[$aa2]}}/$maximum;
    }
}
my $trial=20;
my $iteration=200;
my $B=7;
my $D=2;
my $N=61;
    for(my $i=0; $i<$N; $i++){
	for(my $j=0; $j<$N; $j++){
	    $V{$i}{$j}=0.1*rand(1);
	}
    }


my $T0=0;
    my $T=$T0=0.6;
    while($T>=0.4){
		my $t=0;
		$T=$T0=$T0/(2**0.33);
		for(my $k=0; $k<$iteration; $k++){
		    $t=$t+0.005;
		    for(my $x=0; $x<$N; $x++){
			$sum_B=0;
			for(my $i=0; $i<$N; $i++){
			$Boltz{$x}{$i}=exp(-del_energy(\%d,\%V,$x,$i)/$T);

			$sum_B=$sum_B+$Boltz{$x}{$i};
			}
		    

		    for(my $i=0; $i<$N; $i++){
			if($sum_B!=0){
			$V{$x}{$i}=$Boltz{$x}{$i}/$sum_B;
			}
			}
		    }
		    my @visit;
		    my $len=tour(\%d,\%V, \@visit);
		$T=$T0/(1.0+$t);
		    print STDERR $T."\n";
		}

}
my @visit;
my $len=tour(\%d, \%V, \@visit);
print "Tour Length:".$len*$maximum."\n";
print "Codon  Amino Acid\n";
for(my $i=0; $i<$N; $i++){
    print $i." ".$codons1[$i]." ".$codon_table{$codons1[$visit[$i]]}."\n";
}
sub get_neighbors{
    my $source=shift;
    my ($first,$second,$third)=split //, $source;
    my @neighbors;
    my @raw_bases=qw(A T C G);
    my %transition=("A"=>"G", "G"=>"A", "C"=>"T", "T"=>"C");
    foreach my $rawbase(@raw_bases){
	next if $rawbase eq $third;
	push(@neighbors, $first.$second.$rawbase);
    }
    return \@neighbors;
}
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
    my $len=0;
    my $ord=0;
    for(my $x=0; $x<$N; $x++){
	for(my $i=0; $i<$N; $i++){
	    if($$Vref{$x}{$i}>=0.9){
		$$visitref[$i]=$x;
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

__AA_SCALED__
aa      polarity        hydropathy      volume  iso.electricity
X       -1.5993 1.672   -1.5833 -2.6562
A       -0.0796 0.6402  -0.9414 0.1086
R       0.7182  -1.391  0.9842  2.302
D       2.2     -1.0686 -0.4652 -1.3798
N       1.0602  -1.0686 -0.4238 -0.1633
C       -0.9154 0.8659  -0.4445 -0.3199
E       2.01    -1.0686 -1.4115 -1.1724
Q       0.5283  -1.0686 0.1767  -0.0527
G       0.2623  -0.0691 -1.5212 0.0948
H       0.4523  -0.9719 0.4045  0.8413
I       -0.8774 1.5108  0.7151  0.1178
L       -0.8774 1.2851  0.7151  0.0994
K       1.0982  -0.9074 0.8807  1.832
M       -0.7255 0.6725  0.5908  -0.0112
F       -0.8395 0.9627  1.1499  -0.131
P       -0.2316 -0.456  -0.9104 0.2469
S       0.1104  -0.1981 -0.9207 -0.0388
T       -0.2316 -0.1658 -0.3203 0.1823
W       -0.7635 -0.2303 1.9367  0.0579
Y       -0.6875 -0.3593 1.2327  -0.0481
V       -0.6115 1.4141  0.156   0.0902

__CODON_FILE__
TTT
TTC
TTA
TTG
CTG
CTA
CTT
CTC
ATC
ATT
ATA
ATG
GTG
GTA
GTT
GTC
GCC
GCT
GCA
GCG
ACG
ACA
ACT
ACC
TCC
TCT
TCA
TCG
CCG
CCA
CCT
CCC
CAC
CAT
CAA
CAG
TAG
TAA
TAT
TAC
AAC
AAT
AAA
AAG
GAG
GAA
GAT
GAC
GGC
GGT
GGA
GGG
AGG
AGA
AGT
AGC
TGC
TGT
TGA
TGG
CGG
CGA
CGT
CGC

__AA_AS_CITIES__
T
G
R
I
R
R
G
Y
T
K
V
K
Y
L
F
P
L
W
L
Q
X
R
X
D
S
P
A
C
H
P
S
S
I
A
R
E
S
G
G
S
L
V
C
X
P
F
Q
V
T
E
R
H
A
M
T
L
N
V
L
D
I
S
A
N
