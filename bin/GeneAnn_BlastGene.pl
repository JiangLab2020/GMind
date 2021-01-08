###   此程序用于根据蛋白质序列搜索基因组中的同源基因

use strict;
use warnings;
use Bio::Perl;

my $InputFas = $ARGV[0];     ###    输入基因组Fasta序列
my $InputPep = $ARGV[1];     ###    输入需要查询的蛋白序列

my @b=split(/[\/\.]/,$InputFas); my $Target=$b[$#b-1];
my $Gene_Modeller=$ENV{'Gene_Modeller'};
my $Gene_Modeller_lib=$ENV{'Gene_Modeller_lib'};

my $OutFlag = 0;
if ( $ARGV[2] ) {
	$OutFlag = $ARGV[2];
}

my $CPU_Num = 20;
if ( $ARGV[3] ) {
	$CPU_Num = $ARGV[3];
}

my $BlastNum = 1;

my $MaxGeneLen = 500;

my $MaxUtrLen  = 500;
if ( $ARGV[4] ) {
	$MaxUtrLen = $ARGV[4];
}

my $MaxIntronLen = 200;       #  内含子最大间隔

if ( !(-e "$OutFlag\.txt") ) {
	`perl $Gene_Modeller/Blast_MultiThreads.pl $InputPep $InputFas $OutFlag\.txt 3 $BlastNum $CPU_Num`;
}

`perl $Gene_Modeller/Blast_Uni.pl $OutFlag\.txt`;

my @line=();
open INPUT,"$OutFlag\_Uni.txt";
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_); 
	my @b = (); 
		
	$b[0] = $a[1]; $b[1] = $a[8]; $b[2] = $a[9]; $b[3] = "+"; 
	$b[4] = $a[0]; $b[5] = $a[6]; $b[6] = $a[7]; $b[7] = $a[2]; $b[8] = $a[3]; $b[9] = $a[10];

	if ( $a[6] > $a[7] ) {
		$b[5] = $a[7]; $b[6] = $a[6]; $b[3] = "-";
	}

	push @line,[@b];
}
close(INPUT);

my %QList = ();
my %DList = ();

my %AList = ();  my $Num = 0;

open OUTPUT,">$OutFlag\.sort";

if ( @line ) {	
	my @linesort=sort { $a->[4] cmp $b->[4] or $a->[0] cmp $b->[0] or $a->[5] <=> $b->[5] } @line;
	foreach ( @linesort ) {		
		my @a = @$_;
		push @{$QList{$a[4]."\t".$a[0]."\t".$a[3]}},$a[1]."_".$a[2];
		push @{$DList{$a[4]."\t".$a[0]."\t".$a[3]}},$a[5]."_".$a[6];

		print OUTPUT join("\t",@$_),"\n";
	}
}

close(OUTPUT);

my %CList = ();  my %GList = ();

open OUTPUT1,">$OutFlag\.sort.link";
open OUTPUT2,">$OutFlag\.sort.rem";
open OUTPUT3,">$OutFlag\.sort.cut";           ###   记录了哪些区域需要把基因组划分开(UTR),记录哪些区域需要倍Mask(Intron)
open OUTPUT4,">$OutFlag\.sort.partial";
open OUTPUT5,">$OutFlag\.sort.cutno";

foreach my $key (sort keys %DList) {
	my @a = @{$DList{$key}};
	my @Gene = @{$QList{$key}};
	my @Name = split(/\t/,$key);

	my @d = (); my $k = 0;  my $Fomer = 0;  my @f = ();  my %GapList = ();

	my $MinGene = $MaxGeneLen; my $MaxGene = 0;

	for (my $i=0;$i<@Gene ;$i++) {
		my @e = split(/\_/,$Gene[$i]);
		if ( $e[0] < $MinGene ) {
			$MinGene = $e[0];
		}
		if ( $e[1] > $MaxGene ) {
			$MaxGene = $e[1];
		}
		if ( $i > 0 ) {
			if ( $Name[2] eq "+" and $e[0] < $Fomer ) {
				$k++;
				$GapList{$i-1}+=1;
			}
			if ( $Name[2] eq "-" and $e[0] > $Fomer ) {
				$k++;
				$GapList{$i-1}+=1;
			}
		}

		push @{$d[$k]},$Gene[$i];
		push @{$f[$k]},$a[$i];

		$Fomer = $e[0];
	}

	if ( ($MaxGene-$MinGene) > $MaxGeneLen*0.6 ) {
		my @GeneLim = ();
		for (my $i=0;$i<$k+1 ;$i++) {
			print OUTPUT1 $key,"\t",join("\t",@{$f[$i]}),"\t",join("\t",@{$d[$i]}),"\t",$i+1,"\n";

			$MinGene = $MaxGeneLen;   $MaxGene = 0;
			my @b = @{$d[$i]};
			for (my $j=0;$j<@b ;$j++) {
				my @e = split(/\_/,$b[$j]);
				if ( $e[0] < $MinGene ) {
					$MinGene = $e[0];
				}
				if ( $e[1] > $MaxGene ) {
					$MaxGene = $e[1];
				}
			}
			push @GeneLim,"$MinGene\_$MaxGene";
		}

		$GList{$key}=join("\t",@GeneLim);

		###################################################################################################      定义哪段区域位为UTR, 哪段区域为内含子
		my $EndLast = 0;
		if ( $Name[0] =~ /.*\_(\d+)\_(\d+)/ ) {
			$EndLast = $2-$1+1;
		}

		my @LocRem = ();  my @Loc = @a;  my @GeneRem = ();  my @GapMark = ();
		for (my $i=0;$i<@Loc ;$i++) {
			my @b = split(/\_/,$Loc[$i]);
			my @bg= split(/\_/,$Gene[$i]);

			if ( $i == 0 and $b[0] > 1 ) {
				my $Beg = 1; my $End = $b[0]-1;
				push @LocRem,"$Beg\_$End";

				   $Beg = 1;    $End = $bg[0]-1;
				if ( $Name[2] eq "-" ) {
					$Beg = $bg[1]+1;   $End = $MaxGeneLen;
				}
				push @GeneRem,"$Beg\_$End";
				push @GapMark,"UTR";
			}
			if ( $i < $#Loc ) {
				my @c = split(/\_/,$Loc[$i+1]);
				my $Beg = $b[1]+1; my $End = $c[0]-1;
				push @LocRem,"$Beg\_$End";

				my @cg= split(/\_/,$Gene[$i+1]);
				   $Beg = $bg[1]+1;   $End = $cg[0]-1;
				if ( $Name[2] eq "-" ) {
					$Beg = $cg[1]+1;    $End = $bg[0]-1;
				}
				push @GeneRem,"$Beg\_$End";
				if ( $GapList{$i} ) {
					push @GapMark,"UTR";
				}
				else {
					push @GapMark,"Intron";
				}
			}
			if ( $i == $#Loc and $b[1] < $EndLast ) {
				my $Beg = $b[1]+1; my $End = $EndLast;
				push @LocRem,"$Beg\_$End";

				   $Beg = $bg[1]+1;   $End = $MaxGeneLen;
				if ( $Name[2] eq "-" ) {
					$Beg = 1;    $End = $bg[0]-1;
				}
				push @GeneRem,"$Beg\_$End";
				push @GapMark,"UTR";
			}
		}

		print OUTPUT2 $key,"\t",join("\t",@Loc),"\tXXX\t",join("\t",@Gene),"\n";

		my @MarkLocRem = ();   my @MarkGeneRem = ();  my @Cut = ();

		for (my $i=0;$i<@LocRem ;$i++) {
			push @MarkLocRem,$GapMark[$i]."|".$LocRem[$i];
			push @MarkGeneRem,$GapMark[$i]."|".$GeneRem[$i];

			my @b = split(/\_/,$LocRem[$i]);
			my @c = split(/\_/,$GeneRem[$i]);
			my $LenLoc  = $b[1]-$b[0]+1;
			my $LenGene = $c[1]-$c[0]+1;

			if ( $GapMark[$i] eq "UTR" ) {
				if ( ($i == 0 and $#LocRem > 0) or ($#LocRem == 0 and $b[0] == 1) ) {                ###   第一个UTR或者只有一个左侧的UTR
					if ( $LenLoc > $MaxUtrLen ) {
						my $Beg = $b[0]; my $End = $b[1]-$MaxUtrLen;
						push @Cut,"UTR\|$Beg\_$End";
					}
				}
				elsif ( $i == $#LocRem ) {                                        ###   最后一个UTR
					if ( $LenLoc > $MaxUtrLen ) {
						my $Beg = $b[0]+$MaxUtrLen; my $End = $b[1];
						push @Cut,"UTR\|$Beg\_$End";
					}
				}
				else {
					if ( $LenLoc > $MaxUtrLen*2 ) {
						my $Beg = $b[0]+$MaxUtrLen; my $End = $b[1]-$MaxUtrLen;
						push @Cut,"UTR\|$Beg\_$End";
					}
					elsif ( $LenLoc > 1 ) {
						my $Beg = $b[0]+int($LenLoc/2); my $End = $b[1]-int($LenLoc/2);
						push @Cut,"UTR\|$Beg\_$End";
					}
				}
			}
			if ( $GapMark[$i] eq "Intron" ) {				
				if ( $LenLoc > $MaxIntronLen and $LenGene < 5 ) {           ###   400代表基因组上内含子长度,  5代表两个内含子中可能包含的外显子长度
					my $Beg = $b[0]+int($MaxIntronLen/2); my $End = $b[1]-int($MaxIntronLen/2);
					push @Cut,"Intron\|$Beg\_$End";
				}
			}
		}
		if ( @Cut ) {
			print OUTPUT3 $key,"\t",join("\t",@Cut),"\n";
			$CList{$key}=join("\t",@Cut);
		}
		else {
			print OUTPUT5 $key,"\t",join("\t",@Loc),"\tXXX\t",join("\t",@Gene),"\n";
		}
		print OUTPUT2 $key,"\t",join("\t",@MarkLocRem),"\tXXX\t",join("\t",@MarkGeneRem),"\n";
	}
	else {
		print OUTPUT4 $key,"\t",join("\t",@a),"\tXXX\t",join("\t",@Gene),"\n";
	}

}

close(OUTPUT1);
close(OUTPUT2);
close(OUTPUT3);
close(OUTPUT4);
close(OUTPUT5);

#############################################################################  2步,首先将Intron区域Mask掉,再基于UTR区域将基因组划分开

my %SList = ();
my $seqdata = $InputFas;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) { 
	$SList{$seq->display_id}=$seq->seq();
}

my %LList = ();
$seqdata = $InputPep;
$seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) { 
	$LList{$seq->display_id} = length($seq->seq());
}


open OUTPUT1,">$Target\_Cut_$MaxUtrLen\_$MaxIntronLen\.fasta";
open OUTPUT2,">$Target\_Cut_$MaxUtrLen\_$MaxIntronLen\.txt";

foreach my $key (sort keys %CList) {
	my @Name = split(/\t/,$key);
	my @Cut  = split(/\t/,$CList{$key});
	my $MySeq= $SList{$Name[0]};

	my @UTR = ();  my @Intron = ();

	for (my $i=0;$i<@Cut ;$i++) {
		if ( $Cut[$i] =~ /UTR\|(.*)/ ) {
			push @UTR,$1;
		}
		if ( $Cut[$i] =~ /Intron\|(.*)/ ) {
			push @Intron,$1;
		}
	}

	if ( @Intron ) {
		$MySeq = TranSeqN($MySeq,[@Intron]);
	}
	
	my @Loc = @UTR;  my @LocRem = ();

	for (my $i=0;$i<@Loc ;$i++) {
		my @b = split(/\_/,$Loc[$i]);
		if ( $i == 0 and $b[0] > 1 ) {
			my $Beg = 1; my $End = $b[0]-1;
			push @LocRem,"$Beg\_$End";
		}
		if ( $i < $#Loc ) {
			my @c = split(/\_/,$Loc[$i+1]);
			my $Beg = $b[1]+1; my $End = $c[0]-1;
			if ( $Beg < $End ) {
				push @LocRem,"$Beg\_$End";
			}
		}
		if ( $i == $#Loc and $b[1] < length($MySeq) ) {
			my $Beg = $b[1]+1; my $End = length($MySeq);
			push @LocRem,"$Beg\_$End";
		}
	}

	my @New = ();   my @GeneLim = split(/\t/,$GList{$key});
	
	for (my $i=0;$i<@LocRem ;$i++) {
		my @b = split(/\_/,$LocRem[$i]);
		my $NameNew = $Name[0];
		if ( $NameNew =~ /(.*)\_(\d+)\_(\d+)$/) {
			my $Chr = $1;
			my $Beg = $2;
			my $End = $3;
			my $x0 = $b[0];
			my $y0 = $b[1];
			my $x1 = $Beg+$x0-1;
			my $y1 = $x1+$y0-$x0;
			$NameNew = $Chr."_".$x1."_".$y1;
		}
		my $MySeq = substr($MySeq,$b[0]-1,$b[1]-$b[0]+1);
		print OUTPUT1 ">",$NameNew," ",$Name[0],"\n";
		print OUTPUT1 $MySeq,"\n";

		push @New,$NameNew;

		if ( $GeneLim[$i] ) {
			my @Lim = split(/\_/,$GeneLim[$i]);
			print OUTPUT2 $key,"\t",$i+1,"\t",$LocRem[$i],"\t",$NameNew,"\t",$GeneLim[$i],"\t",Float(($Lim[1]-$Lim[0]+1)/$LList{$Name[1]}),"\n";
		}
		else {
			print $key,"\t",$LocRem[$i],"\t",$NameNew,"\n";
		}
	}	
}

close OUTPUT1;
close OUTPUT2;

########################################################################################################

sub SortArray{
	my $Str = shift;
	my @a = @{$Str};
	my @line=();
	for (my $i=0;$i<@a ;$i++) {
		my @b = split(/\_/,$a[$i]);
		push @line,[@b];
	}

	my @c = ();
	if ( @line ) {	
		my @linesort=sort { $a->[0] <=> $b->[0] } @line;
		foreach ( @linesort ) {
			push @c,join("\_",@$_);
		}
	}
	return [@c];
}

sub TranSeqN{
	my ($Str, $LocA) = @_;
	my @Loc = @{$LocA};
	my @a = split(//,$Str);
	for (my $i=0;$i<@Loc ;$i++) {
		my @b = split(/\_/,$Loc[$i]);
		for (my $j=$b[0];$j<=$b[1] ;$j++) {
			$a[$j-1] = "N";
		}
	}
	return join("",@a);
}

############    保留小数点后两位
sub Float{
	my $str=shift;
	if ($str =~ /\./) {
		$str = sprintf "%0.3f",$str;
	}
	return $str;
}
