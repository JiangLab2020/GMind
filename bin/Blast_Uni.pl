###     此程序对Blast的所有正链结果进行归一化分为两步第一步相同区域查询范围进行合并第二步将级相邻区域进行链接

use strict;
use warnings;
use Bio::Perl;

my $BlastFile = $ARGV[0];

my @b=split(/[\/\.]/,$BlastFile); my $Target=$b[$#b-1];

my $Threshold = 0.8;   

my %DList = ();

open INPUT,$BlastFile;
open OUTPUT1,">GeneLocPS.txt";
open OUTPUT2,">GeneLocNS.txt";

while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_); 
	if ( $a[6] < $a[7] ) {
		print OUTPUT1 $a[0],"\t",$a[1],"\t",$a[6],"\t",$a[7],"\t",$a[8],"\t",$a[9],"\t",$a[11],"\n";
		$DList{$a[0]."\t".$a[1]."\t".$a[6]."\t".$a[7]."\t".$a[8]."\t".$a[9]."\t".$a[11]} = $_;
	}
	else {
		print OUTPUT2 $a[0],"\t",$a[1],"\t",$a[7],"\t",$a[6],"\t",$a[8],"\t",$a[9],"\t",$a[11],"\n";
		$DList{$a[0]."\t".$a[1]."\t".$a[7]."\t".$a[6]."\t".$a[8]."\t".$a[9]."\t".$a[11]} = $_;
	}
}

close(INPUT);
close(OUTPUT1);
close(OUTPUT2);

my $Gene_Modeller=$ENV{'Gene_Modeller'};
my $Gene_Modeller_lib=$ENV{'Gene_Modeller_lib'};
`perl $Gene_Modeller/Tab_Rank.pl GeneLocPS.txt GeneLocPS_Sort.txt 1 3 4`;
`perl $Gene_Modeller/Tab_Rank.pl GeneLocNS.txt GeneLocNS_Sort.txt 1 3 4`;

`perl $Gene_Modeller/LocRegionUni.pl GeneLocPS_Sort.txt $Threshold`;
`perl $Gene_Modeller/LocRegionUni.pl GeneLocNS_Sort.txt $Threshold`;

open OUTPUT,">$Target\_Uni.txt";

open INPUT,"GeneLocPS_Sort_uni.txt";
while (<INPUT>) {
	chomp($_);
	print OUTPUT $DList{$_},"\n";
}
close(INPUT);

open INPUT,"GeneLocNS_Sort_uni.txt";
while (<INPUT>) {
	chomp($_);
	print OUTPUT $DList{$_},"\n";
}
close(INPUT);

close(OUTPUT);

#my $QueryGapMinLen = 10;  my $DBGapMinLen = 5;
#
#my @AllLine =();
#open INPUT,"GeneLoc_Sort_uni.txt";
#
#while (<INPUT>) {
#	chomp($_);
#	my @a=split(/\t/,$_); 
#	if ( !@AllLine ) {
#		push @AllLine,$_;
#	}
#	else {
#		my @b = split(/\t/,$AllLine[$#AllLine]);
#		if ( $a[0] eq $b[0] and $a[2]-$b[3] < $QueryGapMinLen and $a[4]-$b[5] < $DBGapMinLen ) {
#			$b[3] = $a[3];
#			$b[5] = $a[5];
#			$AllLine[$#AllLine] = join("\t",@b);
#		}
#		else {
#			push @AllLine,$_;
#		}
#	}
#}
#close(INPUT);
#
#
#open OUTPUT,">$Target\_Uni.txt";
#
#for (my $i=0;$i<@AllLine ;$i++) {
#	print OUTPUT $AllLine[$i],"\n";
#}
#
#close(OUTPUT);

unlink("GeneLocPS.txt");
unlink("GeneLocPS_Sort.txt");
unlink("GeneLocPS_Sort_uni.txt");
unlink("GeneLocPS_Sort_class.txt");

unlink("GeneLocNS.txt");
unlink("GeneLocNS_Sort.txt");
unlink("GeneLocNS_Sort_uni.txt");
unlink("GeneLocNS_Sort_class.txt");
