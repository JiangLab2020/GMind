use strict;
use warnings;
use Bio::SeqIO;
use List::Util qw/max min/;

my $InputFA =  $ARGV[0];
my $Group   =  $ARGV[1];  

#my $RefSpe  =  "rice";
my $RefSpe  =  "PUGT";
if ( $ARGV[2] ) {
	$RefSpe  = $ARGV[2];
}

#################         coprinus                                 | Coprinus cinereus              
#################         saccharomyces_cerevisiae_S288C           | Saccharomyces cerevisiae       
#################         aspergillus_fumigatus                    | Aspergillus fumigatus          

#################         human                                    | Homo sapiens                   
#################         fly                                      | Drosophila melanogaster             


my @b=split(/[\/\.]/,$InputFA); my $species=$b[$#b-1];
my @d=split(/\_/,$InputFA);
my $specie_num=$d[0];

my @ref = split(//,$RefSpe);   my $refTar = $ref[0].$ref[1];
my $Abbr="";
my @sp=split(/\_/,$species);
$Abbr=$sp[0]."_aug".$refTar;

my %NList=();	
print "####  1 Counting the number of sequences\n";
########################################################################################################
my $Len=0;
my %AList = ();            
my $specie_name="";				
my $seqdata = $InputFA;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {   
	$AList{$seq->display_id}=$seq->seq();
	$Len++;
}
print "       ","Len = ","$Len","\n";
########################################################################################################

print "####  2 Spliting the sequences to $Group groups\n";
########################################################################################################
my @GroLim=();  my @Cut=();
my $GroLen=int $Len/$Group;
if ($GroLen<1) {
	$Group=$Len;
	$GroLen=1;
}
for (my $i=0;$i<$Group-1 ;$i++) {
	$GroLim[$i]=($i+1)*$GroLen;
	$Cut[$i]=$GroLen;
}
$GroLim[$Group-1]=$Len;
$Cut[$Group-1]=$GroLim[$Group-1]-$GroLim[$Group-2];

for (my $i=0;$i<@Cut ;$i++) {
	print "       ",$i+1,"\t",$Cut[$i],"\t",$GroLim[$i],"\n";
}

my $k=0; my $kk=0; my $Flag=0;
open INPUT,$InputFA;
while (<INPUT>) {
	if ( $_ =~ /^\>/ and $k < $GroLim[$kk] and $Flag == 0 ) {
		if ( $k > 0 ) {
			close(OUTPUT);
		}
		$Flag=1;		
		open OUTPUT,">$species\_$kk\.fas";
	}

	print OUTPUT $_;
	
	if ($_ =~ /^\>/) {
		$k++;
	}
	
	if ( $_ =~ /^\>/ and $k == $GroLim[$kk] ) {	
		$kk++;
		$Flag=0;
	}
}
close(OUTPUT);	
close(INPUT);
########################################################################################################

print "####  3 MultiThreads Augustus\n";
########################################################################################################
my @childs=();
my $count=0;

for (my $j=0;$j<$Group ;$j++) {
	$count++;                      
	my $pid = fork();			
	if ($pid) {                      
		sleep 1;
		push(@childs, $pid);
	}
	elsif ($pid == 0) {
		print "       ","fork process $count\n";	
		
		`augustus --species=$RefSpe $species\_$j\.fas >$species\_$j\.gff`;

		exit(0);
	}
	else {
		die "couldn't fork: $!\n";
	}
}
foreach (@childs) {
	waitpid($_,0);
}
########################################################################################################

print "####  4 Fetching protein sequences from Gff file\n";
########################################################################################################
my $MinLen=100;
my $GeneNum=1;
my $Str="1000000";	

open OUTPUT1,">$species\_aug$refTar\.pep";
open OUTPUT2,">$species\_aug$refTar\.cds";
open OUTPUT3,">$species\_aug$refTar\.gff";
for (my $j=0;$j<$Group ;$j++) {
	my $str="";
	my $strCDS="";
	my $strStrand="";
	my $StartStop="";
	my $ChrName="";
	my $gff_gene="";
	my @gff_cds=();
	my $CDS_start=0;
	my $CDS_end=0;
	my @CDS=();
	open INPUT,"$species\_$j\.gff";
	while (<INPUT>) {
		chomp($_);
		$_=~s/^\# //;
		my @a=split(/\t/,$_);
		if ($a[2] and $a[2] eq "gene" ) {
			$ChrName=$a[0];
			pop @a;
			$gff_gene=join("\t",@a);
		}
		if ( $_ =~ /start gene/ ) {
			$str="";
			$strCDS="";
			$strStrand="";
			$StartStop="";
			$gff_gene="";
			@gff_cds=();
			$CDS_start=0;
			$CDS_end=0;
			@CDS=();
		}
		$str.=$_;
		if ( $_ =~ /end gene .*/ ) {
			if ( $str=~ /\[(.*)\]/) {
				my $SeqPep = $1;
				my $SeqLen = length($1);
				if ( $SeqLen > $MinLen ) {					
					my $NameID=$Abbr."\_".substr($Str,0,length($Str)-length($GeneNum)).$GeneNum; $GeneNum++;
					print OUTPUT3 $gff_gene,"\t","ID=",$NameID,"\n",join("\tID\=$NameID\n",@gff_cds),"\t","ID=",$NameID,"\n";
					$CDS_start=min @CDS; $CDS_end=max @CDS;
					if ( $StartStop =~ /Stop/ ) {
						$SeqPep.="*";
					}
					if ( $NList{$specie_num} ) {
						print OUTPUT1 "\>",$NameID," ",$ChrName," ",$CDS_start,"_",$CDS_end," ",$strStrand," ","[",$NList{$specie_num},"]","\n",$SeqPep,"\n";
					}
					else {
						print OUTPUT1 "\>",$NameID," ",$ChrName," ",$CDS_start,"_",$CDS_end," ",$strStrand," ","[",$species,"]","\n",$SeqPep,"\n";
					}
					
					my $SeqCDS=$strCDS;
					if ( $strStrand =~ /\-/ ) {
						my $seqobj = Bio::Seq->new(-display_id => 'my_id',-seq => $strCDS);
						$SeqCDS=$seqobj->revcom->seq();
					}
					my $strMod = length($SeqCDS) % 3;
					if ( $StartStop =~ /Start/ ) {
						if ( $strMod > 0 ) {
							$SeqCDS=substr($SeqCDS,0,length($SeqCDS)-$strMod);
						}
					}
					elsif ( $StartStop =~ /Stop/ ) {
						if ( $strMod > 0 ) {
							$SeqCDS=substr($SeqCDS,$strMod,length($SeqCDS)-$strMod);
						}						
					}
					else {
						my $SeqPep_10=substr($SeqPep,0,10);
						for (my $i=0;$i<=2 ;$i++) {
							my $LenRem = length($SeqCDS)-$i; 
							my $strMod2= $LenRem % 3;
							my $SeqCDS2=substr($SeqCDS,$i,$LenRem-$strMod2);
							my $SeqPep2=DNA2Pep($SeqCDS2);
							my $SeqPep2_10=substr($SeqPep2,0,10);
							if ( $SeqPep_10 eq $SeqPep2_10 ) {
								$SeqCDS=$SeqCDS2; last;
							}
						}
					}
					if ( $NList{$specie_num} ) {
						print OUTPUT2 "\>",$NameID," ",$ChrName," ",$CDS_start,"_",$CDS_end," ",$strStrand," ","[",$NList{$specie_num},"]","\n",$SeqCDS,"\n";
					}
					else {
						print OUTPUT2 "\>",$NameID," ",$ChrName," ",$CDS_start,"_",$CDS_end," ",$strStrand," ","[",$species,"]","\n",$SeqCDS,"\n";
					}
				}
			}
		}
		if ( $a[2] and $a[2] eq "CDS" ) {
			$strCDS.=FetchLocSeq($AList{$a[0]},$a[3],$a[4],"+");
			$strStrand=$a[6];
			push(@CDS,($a[3],$a[4]));
			pop @a;
			my $gff_c=join("\t",@a);
			push(@gff_cds,$gff_c);
			
		}
		if ( $a[2] and $a[2] eq "start_codon" ) {
			$StartStop.=" Start";
		}
		if ( $a[2] and $a[2] eq "stop_codon" ) {
			$StartStop.=" Stop";
		}
	}
	close(INPUT);
}

close(OUTPUT1);
close(OUTPUT2);
close(OUTPUT3);
for (my $t=0;$t<$Group ;$t++) {
	`rm $species\_$t\.fas`;
	`rm $species\_$t\.gff`;
}
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

sub FetchLocSeq{
	my ($MySeq,$a1,$a2,$Strand) = @_;
	my $str1="";
	if ( $a1 > 0 ) {
		$str1=substr($MySeq,$a1-1,$a2-$a1+1);
		if ( $Strand eq "-" ) {
			my $seqobj = Bio::Seq->new( -display_id => 'my_id',-seq => $str1);
			$str1=$seqobj->revcom->seq();
		}
	}
	return uc($str1);
}


sub DNA2Pep{
	my $DNA = shift;
	my $seqobj = Bio::Seq->new(-display_id => 'my_id',-seq => $DNA);
	my $Pep=$seqobj->translate()->seq();
	return $Pep;
}
