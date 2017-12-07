#!/usr/bin/env perl
# look for optimized genetic code using simulated anneal
# Nov 18, 2017
# Modified from Oliver's AANeighborScaled in following ways
## 1. use aa scaled excluding "X"
## 2. Make "X" very similar to "W"
## 3. Linear codon starting from GTC (highest corr)
## 4. Use -1 & +1, instead of 3 mutational nabes

## Preliminaries
use strict;
use Inline::Files; # necessary to embed multiple data sets
use Data::Dumper;
use Getopt::Std;
use Carp;
my %Boltz;
my $sum_B;
my %V; # activities of individual neurons
my @visit; # house an evolved code
my %opts;
getopts('i:t:b:d:s:', \%opts);
#my $titv = $opts{i} || 4; # transition:transversion ratio; default 4
my %transitions = ('A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C');
my @bases = qw(A T C G);
my %aa_scaled; # aa indices scaled
my %dist_aa; # aa scaled distances
my (@codons, %codon_order);
my @aa_as_cities = qw(A A A A C C D D E E F F G G G G H H I I I K K L L L L L L M N N P P P P Q Q R R R R R R S S S S S S T T T T V V V V W X X X Y Y); # according to Standard Genetic Code (SGC)

my $order = 0;
while(<CODON_LINEAR>){
    chomp;
    push @codons, $_;
    $codon_order{$_} = $order++
}

#print Dumper(\%codon_order); exit;

my @amino_acids;
while(<AA_SCALED>){
    chomp;
    my @aa_data=split;
    next if /^aa/;
    my $aa = $aa_data[0];
    $aa_scaled{$aa}= {
        'pol' => $aa_data[1],
        'hydro' => $aa_data[2],
        'vol' => $aa_data[3],
        'iso' => $aa_data[4]
    };
    push @amino_acids, $aa;
}

#print Dumper(\@codons, \%aa_scaled); exit;

# get aa pairwise distances, include 0's for self pairs
my $max = 0;
for(my $i = 0; $i <= $#amino_acids; $i++){
    my $aa1 = $amino_acids[$i];
    for(my $j = 0; $j <= $#amino_acids; $j++){
        my $aa2 = $amino_acids[$j];
        my $d = sqrt(
        ($aa_scaled{$aa1}->{pol} - $aa_scaled{$aa2}->{pol})**2 +
        ($aa_scaled{$aa1}->{hydro} - $aa_scaled{$aa2}->{hydro})**2 +
        ($aa_scaled{$aa1}->{vol} - $aa_scaled{$aa2}->{vol})**2);
        $dist_aa{$aa1}{$aa2} = $d; # || 0.05; # assign the same aa the smallest distances        
#        warn $d, "\n";
        $max = $d if $d > $max;
    }
}

#warn Dumper(\%dist_aa); exit;
#warn $min, "\t", $max, "\n";
# why normalize? this seems to be the key to generate a valid tour (has to do with D, which is set to be 2 * dmax)
 for(my $i=0; $i<$#amino_acids; $i++){
    my $aa1 = $amino_acids[$i];
    for(my $j=$i+1; $j<$#amino_acids; $j++){
        my $aa2 = $amino_acids[$j];
            $dist_aa{$aa1}{$aa2} = $dist_aa{$aa2}{$aa1} /= $max;     
    }
 }

#print Dumper(\%dist_aa); exit;

#  start simulated annealing
#my $trial = 20;
my $iteration=$opts{i} || 10;
my $B = $opts{b} || 7; # effect ?
my $D = $opts{d} || 2; # effect ?
my $N=64;

# initialize random activities
for(my $i=0; $i<$N; $i++){ # for each aa
	for(my $j=0; $j<$N; $j++){ # for each codon
	    $V{$i}{$j}=0.1*rand(); # why /10?
	}
}

my $T0;
my $T = $T0 = $opts{s} || 0.6; # starting temperature
my $round = 1;
my $len = 0;
my $minT = $opts{t} || 0.4;
while($T >= $minT){ # for each temperature, why 0.4?
    warn "=" x 10, $round++, "\t", $T, "=" x 10, "\n";
    my $t = 0;
    $T = $T0 = $T0/(2**0.33); # why?
    for(my $k = 0; $k < $iteration; $k++){ # for each iteration
        warn "iteration\t", $k, "\n"; 
        $t += 0.005;
        for(my $x = 0; $x < $N; $x++){ # for each aa
            $sum_B=0;
            for(my $i = 0; $ i < $N; $i++){ # for codon
			     $Boltz{$x}{$i} = exp(-&del_energy(\%V, $x, $i)/$T);
			     $sum_B += $Boltz{$x}{$i};
            }

            for(my $i=0; $i<$N; $i++){ # normalize by row
                $V{$x}{$i} = $Boltz{$x}{$i}/$sum_B if $sum_B;
            }
        }
        $len = &tour(\%V, \@visit);
        last if $len == $N; # reached 64
    }
#    &print_tour(\%V);
#    my @visit;
#    my $len=tour(\%d,\%V, \@visit);
#    last if $len;
    $T = $T0/(1.0 + $t);
#    print STDERR $T."\n";
}

for(my $i=0; $i<$N; $i++){
    print $codons[$i], "\t", $aa_as_cities[$visit[$i]], "\n";
}
#print Dumper(\@visit);
# &print_tour(\%V);

exit;

sub get_neighbors{ # return third-base nabes
    my $source=shift;
    my ($first,$second,$third)=split //, $source;
    my @neighbors;
    foreach my $rawbase (@bases){
	   next if $rawbase eq $third;
	   push @neighbors, $first.$second.$rawbase;
    }
    return \@neighbors;
}

sub del_energy{
    my $Vref = shift;
    my $x = shift; # AA(X)
    my $i = shift; # codon index, 0-indexed
    my $term1 = 0;
    my $term2 = 0;
    for(my $y = 0; $y < $N; $y++){
	    $term1 += $$Vref{$y}{$i} unless $x == $y;
    }
    
    my $aa_x = $aa_as_cities[$x];
    
    for(my $y=0; $y<$N; $y++){ # for AA(Y)
        next if $x == $y;
        my $aa_y = $aa_as_cities[$y];
#        next if $aa_x eq $aa_y;
#	    my $neighborref = get_neighbors($codons[$i]);
#	       for(my $k = 0; $k < @{$neighborref}; $k++){
        if ($i == 63) { # right end
            $term2 += $dist_aa{$aa_x}{$aa_y} * ($$Vref{$y}{$i-1} + $$Vref{$y}{0});       
        } elsif ($i == 0) { # left end
            $term2 += $dist_aa{$aa_x}{$aa_y} * ($$Vref{$y}{$N-1} + $$Vref{$y}{$i+1});       
        } else {
            $term2 += $dist_aa{$aa_x}{$aa_y} * ($$Vref{$y}{$i-1} + $$Vref{$y}{$i+1});
	    }
	}    
    return($B * $term1 + $D * $term2);
}

sub tour{
    my $Vref=shift;
    my $visitref=shift;
    my $length = 0;
    my $ord = 0;
    my $warn_str ='';
    my $energy = 0;
    for(my $x = 0; $x < $N; $x++){
	   for(my $i = 0; $i < $N; $i++){
            my $status = $$Vref{$x}{$i} > 0.9 ? 1 : 0;
            $ord += $status;
            $visitref->[$i] = $x if $status; # evolved code
            } 
    }

    $warn_str .= "stops\t" . $ord. "\t";
#    if ($ord != $N){ # not valid tour (has holes)
#	   return 0;
#    }

    my $term1 = 0;
    my $term2 = 0;
    for(my $i=0; $i<$N-1; $i++){
	   $length += $dist_aa{$aa_as_cities[$$visitref[$i]]}{$aa_as_cities[$$visitref[$i+1]]};

       for(my $x=0; $x<$N; $x++){
        my $aa_x = $aa_as_cities[$x];
        for(my $y=0; $y<$N; $y++){ # for AA(Y)
            next if $x == $y;
            my $aa_y = $aa_as_cities[$y];
#            next if $aa_x eq $aa_y;
            $term1 += $$Vref{$y}{$i} * $$Vref{$x}{$i};
         }
        }
    }

    for(my $x=0; $x<$N; $x++){
        my $aa_x = $aa_as_cities[$x];
        for(my $y=0; $y<$N; $y++){ # for AA(Y)
            next if $x == $y;
            my $aa_y = $aa_as_cities[$y];
            next if $aa_x eq $aa_y;
            for(my $i=0; $i<$N-1; $i++) {            
                if ($i == 63) { # right end
                    $term2 += $dist_aa{$aa_x}{$aa_y} * ($$Vref{$y}{$i-1} + $$Vref{$y}{0});       
                } elsif ($i == 0) { # left end
                    $term2 += $dist_aa{$aa_x}{$aa_y} * ($$Vref{$y}{$N-1} + $$Vref{$y}{$i+1});       
                } else {
                    $term2 += $dist_aa{$aa_x}{$aa_y} * ($$Vref{$y}{$i-1} + $$Vref{$y}{$i+1});
	           }
	       }    
       }
    }
    
    $length += $dist_aa{$aa_as_cities[$$visitref[$N-1]]}{$aa_as_cities[$$visitref[0]]};
     
    $energy = $B * $term1 + $D * $term2;
    
    $warn_str .= "tour length = " . sprintf("%.2f", $length) . "\tenergy = " . sprintf("%.2f", $energy);
    warn $warn_str, "\n";
    return $ord;
}

__AA_SCALED__
aa     polarity    hydropathy    volume   iso.electricity
X      -0.9     -0.2  2.0  -0.1
A       -0.1672 0.7638  -1.067  -0.02972
C       -1.043  1.002   -0.5477 -0.5561
D       2.221   -1.039  -0.5694 -1.858
E       2.022   -1.039  -1.559  -1.603
F       -0.9631 1.104   1.12    -0.3241
G       0.191   0.01531 -1.674  -0.0467
H       0.39    -0.9373 0.3402  0.8703
I       -1.003  1.682   0.6651  -0.0184
K       1.067   -0.8693 0.8384  2.087
L       -1.003  1.444   0.6651  -0.04104
M       -0.8437 0.7978  0.5352  -0.1769
N       1.027   -1.039  -0.5261 -0.3637
P       -0.3264 -0.393  -1.035  0.1401
Q       0.4696  -1.039  0.102   -0.2278
R       0.6686  -1.38   0.9466  2.665
S       0.03184 -0.1208 -1.046  -0.2108
T       -0.3264 -0.08676        -0.4178 0.06085
V       -0.7243 1.58    0.08035 -0.05236
W       -0.8835 -0.1548 1.943   -0.09198
Y       -0.8039 -0.2909 1.207   -0.2222
__CODON_LINEAR__
GTC
GTT
GTA
GTG
ATG
ATA
ATT
ATC
TTC
TTT
TTA
TTG
CTG
CTA
CTT
CTC
CCC
CCT
CCA
CCG
TCG
TCA
TCT
TCC
ACC
ACT
ACA
ACG
GCG
GCA
GCT
GCC
GAC
GAT
GAA
GAG
AAG
AAA
AAT
AAC
TAC
TAT
TAA
TAG
CAG
CAA
CAT
CAC
CGC
CGT
CGA
CGG
TGG
TGA
TGT
TGC
AGC
AGT
AGA
AGG
GGG
GGA
GGT
GGC