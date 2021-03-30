use strict;
use warnings;
use Bio::Perl;

my $InputFas = $ARGV[0];     
my $InputPep = $ARGV[1];     

my @b=split(/[\/\.]/,$InputFas); my $Target=$b[$#b-1];

my $OutFlag = 0;
if ( $ARGV[2] ) {
	$OutFlag = $ARGV[2];
}

my $CPU_Num = 20;
if ( $ARGV[3] ) {
	$CPU_Num = $ARGV[3];
}
my $Genome_dir = 0;
if ( $ARGV[4] ) {
	$Genome_dir = $ARGV[4];
}

my $GMind=$ENV{'GMind'};
#my $Identity = 10;

my $CutLen   = 20;

my $BlastNum = 1000;

my $MaxInterval = 20000;   

my $ExtandLen   = 6000;

my $EValue = 1e-6;

my $GeneLen = 300;

my @args = stat ($InputFas);
my $size = $args[7]/1e+6;

print "The genome size is: ",$size,"M","\n";

if ( $size < 100 ) {
	$MaxInterval = 10000;
}
elsif ( $size < 500 ) {
	$MaxInterval = 20000;
}
elsif ( $size < 1000 ) {
	$MaxInterval = 30000;
}
elsif ( $size < 2000 ) {
	$MaxInterval = 40000;
}
else {
	$MaxInterval = 50000;
}

print "The max intron length is smaller than: ",$MaxInterval,"\n";

if ( !(-e "$OutFlag\.txt") ) {
	`perl $GMind/Blast_MultiThreads.pl $InputFas $InputPep $OutFlag\.txt 4 $BlastNum $CPU_Num`;
}

my @line=();
open INPUT,"$OutFlag\.txt";
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_); 
	if ( $a[10] <= $EValue ) {
		my @b = (); 
		
		$b[0] = $a[1]; $b[1] = $a[8]; $b[2] = $a[9]; $b[3] = "+"; 
		$b[4] = $a[0]; $b[5] = $a[6]; $b[6] = $a[7]; $b[7] = $a[2]; $b[8] = $a[3]; $b[9] = $a[10];

		if ( $a[8] > $a[9] ) {
			$b[1] = $a[9]; $b[2] = $a[8]; $b[3] = "-";
		}

		push @line,[@b];
	}
}
close(INPUT);

my %AList = ();  my $Num = 0;

open OUTPUT,">$OutFlag\.sort";

if ( @line ) {	
	my @linesort=sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] } @line;
	foreach ( @linesort ) {		
		my @a = @$_;
		if ( !$AList{$a[0]} ) {
			push @{$AList{$a[0]}},"$a[1]\t$a[2]";
			$Num = 0;
		}
		else {
			my @b = @{$AList{$a[0]}}; my $Flag = 0;
			for (my $i=0;$i<@b ;$i++) {
				if ( RegionCombine($b[$i],"$a[1]\t$a[2]") ) {
					$AList{$a[0]}[$i] = RegionCombine($b[$i],"$a[1]\t$a[2]");
					$Num = $i;
					$Flag = 1;
					last;
				}
			}
			if ( $Flag == 0 ) {
				push @{$AList{$a[0]}},"$a[1]\t$a[2]";
				$Num = $#b+1;
			}
		}
		print OUTPUT $Num,"\t",join("\t",@$_),"\n";
	}
}

close(OUTPUT);


my %BList = ();

open INPUT,"$OutFlag\.sort";
while (<INPUT>) {
	my @a=split(/\t/,$_);
	my $key = $a[0]."\t".$a[1]."\t".$a[5];
	if ( !$BList{$key} ) {
		push @{$BList{$key}},"$a[6]\t$a[7]";
	}
	else {
		my @b = @{$BList{$key}}; my $Flag = 0;
		for (my $i=0;$i<@b ;$i++) {
			if ( RegionCombine($b[$i],"$a[6]\t$a[7]") ) {
				$BList{$key}[$i] = RegionCombine($b[$i],"$a[6]\t$a[7]");
				$Flag = 1;
				last;
			}
		}
		if ( $Flag == 0 ) {
			push @{$BList{$key}},"$a[6]\t$a[7]";
		}
	}
}
close(INPUT);

my %CList = (); my %DList = ();

open OUTPUT,">$OutFlag\.sort.num";

foreach my $key (sort keys %BList) {
	my @b = @{$BList{$key}};
	for (my $i=0;$i<@b ;$i++) {		
		my @c = split(/\t/,$key);
		my @d = split(/\t/,$b[$i]);
		if ( !$CList{$c[0]."\t".$c[1]} or $DList{$c[0]."\t".$c[1]} < $d[1]-$d[0] ) {
			$CList{$c[0]."\t".$c[1]} = $c[2]."\t".$b[$i];
			$DList{$c[0]."\t".$c[1]} = $d[1]-$d[0];
		}
		print OUTPUT $key,"\t",$b[$i],"\n";
	}
}

close(OUTPUT);


my %EList = ();  my $RegionNum = 0;

open OUTPUT,">$OutFlag\.sort.combine";
open OUTPUT1,">$OutFlag\.sort.combine.short";


foreach my $key (sort keys %AList) {
	my @b = @{$AList{$key}};
	for (my $i=0;$i<@b ;$i++) {	
		$RegionNum++;
		my @c = split(/\t/,$CList{$i."\t".$key});
		if ( $c[2]-$c[1] > $GeneLen ) {
			print OUTPUT $i,"\t",$key,"\t",$b[$i],"\t",$CList{$i."\t".$key},"\n";
			$EList{$i."\t".$key} = $c[2]-$c[1];
		}
		else {
			print OUTPUT1 $i,"\t",$key,"\t",$b[$i],"\t",$CList{$i."\t".$key},"\n";
		}
	}
}

close(OUTPUT);
close(OUTPUT1);

print "All homologs regions is: ",$RegionNum,"\n";

open OUTPUT,">$Target\_$OutFlag\.fasta";

my $seqdata = $InputFas;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	if ( $AList{$seq->display_id} ) {
		my @b = @{$AList{$seq->display_id}};
		for (my $i=0;$i<@b ;$i++) {
			if ( $EList{$i."\t".$seq->display_id} ) {
				my @c = split(/\t/,$b[$i]);

				my $GeneBig = $c[0] - $ExtandLen;
				my $GeneEnd = $c[1] + $ExtandLen;

				if ( $EList{$i."\t".$seq->display_id} < 455 ) {
					$GeneBig = $c[0] - 15000;
					$GeneEnd = $c[1] + 15000;
				}

				if ( $GeneBig < 1 ) {
					$GeneBig = 1;
				}				
				if ( $GeneEnd > length($seq->seq()) ) {
					$GeneEnd  = length($seq->seq());
				}

				my $str1 = substr($seq->seq(),$GeneBig-1,$GeneEnd-$GeneBig+1);

				print OUTPUT ">",$seq->display_id."_".$GeneBig."_".$GeneEnd,"\n";
				print OUTPUT $str1,"\n";
			}
		}
	}
}

close(OUTPUT);

if ($Genome_dir) {
	`cp $Target\_$OutFlag\.fasta $Genome_dir`;
}

###############################################################################################################################

########################################################     获得两个区域的重叠区域的长度
sub RegionOverlap{
	my ($Str1,$Str2) = @_;
	my @a = split(/\t/,$Str1);
	my @b = split(/\t/,$Str2);
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

########################################################     获得两个区域的并集
sub RegionCombine{
	my ($Str1,$Str2) = @_;
	my @a = split(/\t/,$Str1);
	my @b = split(/\t/,$Str2);
	if ( $a[1] < $b[0]-$MaxInterval or $b[1] < $a[0]-$MaxInterval ) {
		return 0;
	}
	else {
		my $c0 = $a[0];
		if ( $b[0] < $a[0] ) {
			$c0 = $b[0];
		}
		my $c1 = $a[1];
		if ( $b[1] > $a[1] ) {
			$c1 = $b[1];
		}
		return "$c0\t$c1";
	}
}
