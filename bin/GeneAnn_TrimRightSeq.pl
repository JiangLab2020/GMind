use strict;
use warnings;
use Bio::SeqIO;

my $TargetDir  = $ARGV[0];   $TargetDir=~s/\/$//;    ####  基因组所在文件夹

my $divCut = 0.15;

opendir DIR,"$TargetDir";
my @dots = readdir(DIR);
close(DIR);

chdir $TargetDir;

my $NumRight = 0;  my $NumWrong = 0;

my %TList = ();    ###         记录哪些是正确注释的基因

open IN,"All_uni_DB_adj_div.txt";
while (<IN>) {
	chomp($_);
	my @a = split(/\t/,$_);
	if ( $a[1] < $divCut ) {
		$TList{$a[0]} = 1+$a[1];		
	}
}
close(IN);

open OUTPUT1,">All_uni_DB_adj_div_good.pep";
open OUTPUT2,">All_uni_DB_adj_div_bed.pep";

my %LList = ();    ###         记录每个染色体所有正确基因的起始位置

my $seqdata = "All_uni_DB_adj_div.pep";
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	my @a = split(/ /,$seq->description);
	if ( $TList{$seq->display_id} ) {
		print OUTPUT1 ">",$seq->display_id," ",$seq->description,"\n";
		print OUTPUT1 $seq->seq(),"\n";
		$NumRight++;
		my @b = split(/\_/,$a[1]);
		$b[0]=$b[0]-1;
		$b[1]=$b[1]-1;
		push @{$LList{$a[0]}},"$b[0]\_$b[1]";
	}
	else {
		print OUTPUT2 ">",$seq->display_id," ",$seq->description,"\n";
		print OUTPUT2 $seq->seq(),"\n";
		$NumWrong++;
	}
}

close OUTPUT1;
close OUTPUT2;

print "RightNum = ",$NumRight,"  ","WrongNum = ",$NumWrong,"\n";

open OUTPUT2,">Trim_All.fasta";
open OUTPUT3,">Trim_mask.txt";

for (my $j=2;$j<@dots ;$j++) {
	if ( $dots[$j] =~ /\.fasta/ and $dots[$j] !~ /^Trim/ ) {
		my $seqdata = $dots[$j];
		my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
		while ( my $seq = $seqio->next_seq ) { 
			if ( !$LList{$seq->display_id} ) {
				print OUTPUT2 ">",$seq->display_id,"\n";
				print OUTPUT2 $seq->seq(),"\n";
			}
			else {
				my @Loc = @{SortArray($LList{$seq->display_id})};  print OUTPUT3 $seq->display_id,"\t",join("\t",@Loc),"\n";
				my @LocRem = ();
				for (my $i=0;$i<@Loc ;$i++) {
					my @b = split(/\_/,$Loc[$i]);
					if ( $i == 0 and $b[0] > 0 ) {
						my $Beg = 0; my $End = $b[0]-1;
						push @LocRem,"$Beg\_$End";
					}
					if ( $i < $#Loc ) {
						my @c = split(/\_/,$Loc[$i+1]);
						my $Beg = $b[1]+1; my $End = $c[0]-1;
						if ( $Beg < $End ) {
							push @LocRem,"$Beg\_$End";
						}
						else {                             ###      注释正确基因出现重叠的情况，重新注释
							$Beg = $b[0]; $End = $c[1];
							push @LocRem,"$Beg\_$End";
						}
					}
					if ( $i == $#Loc and $b[1] < length($seq->seq())-1 ) {
						my $Beg = $b[1]+1; my $End = length($seq->seq())-1;
						push @LocRem,"$Beg\_$End";
					}
				}

				for (my $i=0;$i<@LocRem ;$i++) {
					my @b = split(/\_/,$LocRem[$i]);
					my $Name = $seq->display_id;
					if ( $Name =~ /(.*)\_(\d+)\_(\d+)$/) {
						my $Chr = $1;
						my $Beg = $2;
						my $End = $3;
						my $x0 = $b[0]+1;
						my $y0 = $b[1]+1;
						my $x1 = $Beg+$x0-1;
						my $y1 = $x1+$y0-$x0;
						$Name = $Chr."_".$x1."_".$y1;
					}
					my $MySeq = substr($seq->seq(),$b[0],$b[1]-$b[0]+1);

					print OUTPUT2 ">",$Name," ",$seq->display_id,"\n";
					print OUTPUT2 $MySeq,"\n";
				}
			}
		}
	}
}

close OUTPUT2;
close OUTPUT3;


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
