use lib("/Users/oliverattie/Downloads/Switch-2.17/blib/lib/");
use Switch;
my $exe, $trial, $iteration, $flag;
my $start_time, $elapse_time, $new_E, $t, $T0, $T;
my %Boltz;
my $sum_B;
my %d;
my %distance;
my %V;
my %city;
my @visit;
my @tsp_len;
my @tsp_time;
my $sum_len, $sum_time, $min_E, $max_E, $temp;
$trial=20;
$iteration=1000;
my $B=7;
my $D=2;
my $N=50;
srand(10);
my @codons;
my %codon;
$codons[0]="GGG";
$codons[1]="GGA";
$codons[2]="GGT";
$codons[3]="GGC";
$codons[4]="GAG";
$codons[5]="GAA";
$codons[6]="GAT";
$codons[7]="GAC";
$codons[8]="GTG";
$codons[9]="GTA";
$codons[10]="GTT";
$codons[11]="GTC";
$codons[12]="GCG";
$codons[13]="GCA";
$codons[14]="GCT";
$codons[15]="GCC";
$codons[16]="AGG";
$codons[17]="AGA";
$codons[18]="AGC";
$codons[19]="AAG";
$codons[20]="AAA";
$codons[21]="AAT";
$codons[22]="AAC";
$codons[23]="ATG";
$codons[24]="ATA";
$codons[25]="ATT";
$codons[26]="ATC";
$codons[27]="ACG";
$codons[28]="ACA";
$codons[29]="ACT";
$codons[30]="ACC";
$codons[31]="TGG";
$codons[32]="TGA";
$codons[33]="TGT";
$codons[34]="TGC";
$codons[35]="TAG";
$codons[36]="TAA";
$codons[37]="TAT";
$codons[38]="TAC";
$codons[39]="TTG";
$codons[40]="TTA";
$codons[41]="TTT";
$codons[42]="TTC";
$codons[43]="TCG";
$codons[44]="TCA";
$codons[45]="TCT";
$codons[46]="TCC";
$codons[47]="CGG";
$codons[48]="CGA";
$codons[49]="CGT";
for(my $i=0; $i<$N; $i++){
    $codon{$i}{0}=substr($codons[$i],0,1);
    $codon{$i}{1}=substr($codons[$i],1,1);
    $codon{$i}{2}=substr($codons[$i],2,1);
}
#    print "City ".$i." ".$city{$i}{0}." ".$city{$i}{1}."\n";
#    $flag=1;

#	for(my $j=0; $j<$i; $j++){
#	    if($city{$j}{0}==$city{$i}{0}&&$city{$j}{1}==$city{$i}{1}){
#	    $city{$i}{0}=int(rand(1))*(3.0);
#	    $city{$i}{1}=int(rand(1))*(1.6);
#	    break;
#	}
#    }
#}

for(my $i=0; $i<$N; $i++){
    for(my $j=$i+1; $j<$N; $j++){
	$distance{$i}{$j}=sqrt(transv_dist($codon{$i}{0},$codon{$j}{0})*transv_dist($codon{$i}{0},$codon{$j}{0})+transv_dist($codon{$i}{1},$codon{$j}{1})*transv_dist($codon{$i}{1},$codon{$j}{1})+transv_dist($codon{$i}{2},$codon{$i}{2})*transv_dist($codon{$i}{2},$codon{$j}{2}));
	$distance{$j}{$i}=$distance{$i}{$j};
	print "Distances:".$i." ".$j." ".$distance{$i}{$j}."\n";
    }
}
$maximum=0;
for(my $codon1=0; $codon1<$N; $codon1++){
    for(my $codon2=0; $codon2<$N; $codon2++){
	if($distance{$codon1}{$codon2}>$maximum){
	$maximum=$distance{$codon1}{$codon2};
	}
    }
}
for(my $codon1=0; $codon1<$N; $codon1++){
    for(my $codon2=0; $codon2<$N; $codon2++){
	$d{$codon1}{$codon2}=$distance{$codon1}{$codon2}/$maximum;
    }
}
#for($exe=0; $exe<$trial; $exe++){
    srand(1);
    for(my $i=0; $i<$N; $i++){
	for(my $j=0; $j<$N; $j++){
	    $V{$i}{$j}=0.1*rand(1);
	}
    }
    $T=$T0=1.6;
    while($T>=0.0001){
		$t=0;
		$T=$T0=$T0/(2**0.33);
		for(my $k=0; $k<$iteration; $k++){
		    $t=$t+0.001;
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
#			print "V:".$V{$x}{$i}." ";
			
			}
#			print "\n";
		    }
		    my @visit;
		    $len=tour(\%d,\%V, \@visit);
#		    print "Tour Length:".$len."\n";
#		    print "Stop City X-Coord Y-Coord\n";
#		    for(my $i=0; $i<$N; $i++){
#			print $i." ".$visit[$i]." ".$codon{$visit[$i]}{0}.$codon{$visit[$i]}{1}.$codon{$visit[$i]}{2}."\n";
#		    }
		$T=$T0/(1.0+$t);
		    print STDERR $T."\n";
		}

    }
my @visit;
$len=tour(\%d, \%V, \@visit);
print "Tour Length:".$len*$maximum."\n";
print "Stop Codon Nucleotides\n";
for(my $i=0; $i<$N; $i++){
    print $i." ".$visit[$i]." ".$codon{$visit[$i]}{0}.$codon{$visit[$i]}{1}.$codon{$visit[$i]}{2}."\n";
}
#}
sub transv_dist{
    $nuc1=shift;
    $nuc2=shift;
    $dist=0;
    if($nuc1 eq $nuc2){
	$dist=0;
    }
    if((($nuc1 eq "A")&&($nuc2 eq "G"))||(($nuc1 eq "G")&&($nuc2 eq "A"))){
	$dist=1;
    }
    if((($nuc1 eq "C")&&($nuc2 eq "T"))||(($nuc1 eq "T")&&($nuc2 eq "C"))){
	$dist=1;
    }
    if((($nuc1 eq "A")&&($nuc2 eq "C"))||(($nuc1 eq "C")&&($nuc2 eq "A"))){
	$dist=4;
    }
    if((($nuc1 eq "G")&&($nuc2 eq "T"))||(($nuc1 eq "T")&&($nuc2 eq "G"))){
	$dist=4;
    }
    if((($nuc1 eq "A")&&($nuc2 eq "T"))||(($nuc1 eq "T")&&($nuc2 eq "A"))){
	$dist=4;
    }
    if((($nuc1 eq "C")&&($nuc2 eq "G"))||(($nuc1 eq "G")&&($nuc2 eq "C"))){
	$dist=4;
    }
    return $dist;
}
sub energy{
    $dref=shift;
    $Vref=shift;
    my $term1=0;
    my $term2=0;
    for(my $i=0; $i<$N; $i++){
	for(my $x=0; $x<$N; $x++){
	    for(my $y=0; $y<$N; $y++){
		if($y!=$x){
		    $term1=$term1+$$Vref{$x}{$i}*$$Vref{$y}{$i};
		}
	    }
	}
    }
	for(my $x=0; $x<$N; $x++){
	    for(my $y=0; $y<$N; $y++){
		if($y!=$x){
		    for(my $i=0; $i<$N; $i++){
			switch($i){
			    case 0{$term2=$term2+$$dref{$x}{$y}*$$Vref{$x}{$i}*($$Vref{$y}{$i+1}+$$Vref{$y}{$N-1});break;}
			    case ($N-1){$term2=$term2+$$dref{$x}{$y}*$$Vref{$x}{$i}*($$Vref{$y}{0}+$$Vref{$y}{$i-1});break;}
			  else {$term2=$term2+$$dref{$x}{$y}*$$Vref{$X}{$i}*($$Vref{$y}{$i+1}+$$Vref{$y}{$i-1});break;}

			}
}
}
}
}
			return ($B*0.5*$term1+$D*0.5*$term2);
}		    
	
sub del_energy{
    $dref=shift;
    $Vref=shift;
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
	    switch($i){
		case 0{$term2=$term2+$$dref{$x}{$y}*($$Vref{$y}{$i+1}+$$Vref{$y}{$N-1});break;}

		case ($N-1){$term2=$term2+$$dref{$x}{$y}*($$Vref{$y}{0}+$$Vref{$y}{$i-1});break;}

	      else {$term2=$term2+$$dref{$x}{$y}*($$Vref{$y}{$i+1}+$$Vref{$y}{$i-1});break;}

	    }
	}
    }
	    return($B*$term1+$D*$term2);
}
sub tour{
    $dref=shift;
    $Vref=shift;
    $visitref=shift;
    my $len=0;
    my $ord=0;
    for(my $x=0; $x<$N; $x++){
	for($i=0; $i<$N; $i++){
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
