use strict;
use warnings;
use Bio::Perl;

my $AllMethodPep =  "All_Methods.pep";                     
my $AllMethoddiv =  "All_Methods_div.txt";                 
my $UniMethodPep =  "All_uni_DB.pep";                      

my @b=split(/[\/\.]/,$UniMethodPep); my $Target=$b[$#b-1];

my $divCut = 0.15;  $divCut+=1;

my %DList = ();

open INPUT,$AllMethoddiv;
while (<INPUT>) {
	chomp($_);	
	my @a=split(/\t/,$_);  
	$DList{$a[0]} = $a[1]+1;
}
close(INPUT);

my %LList = ();    
my %TList = ();
my %NList = ();

my $seqdata = $UniMethodPep;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	my @a = split(/ /,$seq->description);
	push @{$LList{$a[0]}},$a[1];
	$TList{$seq->display_id} = $DList{$seq->display_id};
	$NList{$a[0].$a[1]} = $seq->display_id;
}


my %RList = ();        

my $k0 = 0;  my $k1 = 0;  my $k2 = 0;

open OUTPUT,">All_uni_DB_adj.pep";
open OUTPUT1,">All_uni_DB_adj.txt";

$seqdata = $AllMethodPep;
$seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	my @a = split(/ /,$seq->description);
	if ( !$TList{$seq->display_id} and $DList{$seq->display_id} and $DList{$seq->display_id} < $divCut ) {
		my @Loc =  @{$LList{$a[0]}};
		my $MaxOL = 0; my $MaxLoc = "";
		for (my $i=0;$i<@Loc ;$i++) {
			my $OL = RegionOverlap($a[1],$Loc[$i]);
			if ( $OL > $MaxOL ) {
				$MaxOL = $OL; $MaxLoc = $Loc[$i];
			}
		}
		if ( $MaxOL > 0 ) {
			my @c = split(/\_/,$MaxLoc);
			if ( $MaxOL/($c[1]-$c[0]+1) > 0.1 ) {
				if ( $DList{$NList{$a[0].$MaxLoc}} > $divCut and !$RList{$NList{$a[0].$MaxLoc}} ) {
					print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
					print OUTPUT $seq->seq(),"\n";
					$RList{$NList{$a[0].$MaxLoc}} = $seq->display_id;
					print OUTPUT1 $seq->display_id,"\t",$DList{$seq->display_id},"\t",$NList{$a[0].$MaxLoc},"\t",$DList{$NList{$a[0].$MaxLoc}},"\n";
					$k2++;
				}
			}
			else {
				print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
				print OUTPUT $seq->seq(),"\n";
				print OUTPUT1 $seq->display_id,"\t",$DList{$seq->display_id},"\t",$NList{$a[0].$MaxLoc},"\t",$DList{$NList{$a[0].$MaxLoc}},"\t","Less","\n";
				$k1++;
			}
		}
		else {
			print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
			print OUTPUT $seq->seq(),"\n";
			print OUTPUT1 $seq->display_id,"\t",$DList{$seq->display_id},"\n";
			$k0++;
		}
	}
}

my $k3 = 0;

$seqdata = $UniMethodPep;
$seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	if ( !$RList{$seq->display_id} ) {
		$RList{$seq->display_id}+=1;
		print OUTPUT ">",$seq->display_id," ",$seq->description,"\n";
		print OUTPUT $seq->seq(),"\n";
	}
	else {
		$k3++;
	}
}

close(OUTPUT);
close(OUTPUT1);

print $k0," ",$k1," ",$k2," ",$k3,"\n";

open OUTPUT,">All_uni_DB_adj_div.txt";

$seqdata = "All_uni_DB_adj.pep";
$seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	if ( $DList{$seq->display_id} ) {
		print OUTPUT $seq->display_id,"\t",Float($DList{$seq->display_id}-1),"\n";
	}
}

close(OUTPUT);

unlink("All_uni_DB_adj.txt");

########################################################     
sub RegionOverlap{
	my ($Str1,$Str2) = @_;
	my @a = split(/\_/,$Str1);
	my @b = split(/\_/,$Str2);
	if ( $a[1] < $b[0] or $b[1] < $a[0] ) {
		return 0;
	}
	else {
		my $c0 = $a[0];
		if ( $b[0] > $a[0] ) {
			$c0 = $b[0];
		}
		my $c1 = $a[1];
		if ( $b[1] < $a[1] ) {
			$c1 = $b[1];
		}
		return $c1-$c0+1;
	}
}

############
sub Float{
	my $str=shift;
	if ($str =~ /\./) {
		$str = sprintf "%0.3f",$str;
	}
	return $str;
}
