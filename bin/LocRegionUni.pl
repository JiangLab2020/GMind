##   此程序用于获取每个Block涉及的区域大小,以及包含的基因数量

use strict;
use warnings;

my $File = $ARGV[0];
my @b=split(/[\/\.]/,$File); my $Target=$b[$#b-1];

###################################################################################    按照Block重叠区域聚类
my $Threshold = 0.8;
if ( $ARGV[1] ) {
	$Threshold = $ARGV[1];
}

my @Event = ();
my @EventClass = ();

open INPUT,$File;
while (<INPUT>) {
	chomp($_);	
	my @a=split(/\t/,$_); my $Class = 0;   

	my $strA1 = "$a[2]\t$a[3]"; 
	my $LenA1 = $a[3]-$a[2]+1;
	my $ScoreA = $a[$#a];
	my $Flag = 0;
	if ( @Event ) {
		for (my $i=$#Event;$i>$#Event-5 ;$i--) {
			if ( $i >= 0 ) {
				my @c = split(/\t/,$Event[$i]);
				if ( $a[1] eq $c[1] ) {
					my $strC1 = "$c[2]\t$c[3]";
					my $LenC1 = $c[3]-$c[2]+1;
					my $ScoreC = $c[$#c];

					my $OverlapAC1 = RegionOverlap($strA1,$strC1);

					####   两对重叠区域同时大于阈值,则归为一类
					if ( $OverlapAC1/$LenA1 > $Threshold or $OverlapAC1/$LenC1 > $Threshold ) {
						$Flag++;
						if ( $ScoreA > $ScoreC ) {
							$Event[$i]=$_;
						}
						$Class=$i;
						last;
					}
				}
			}
		}
	}
	if ( !@Event or $Flag == 0 ) {
		push @Event,$_;
		$Class=$#Event;
	}
	push @EventClass,$Class."\t".$_;
}
close(INPUT);

open OUTPUT,">$Target\_uni.txt";

for (my $i=0;$i<@Event ;$i++) {
	print OUTPUT $Event[$i],"\n";
#	print OUTPUT $i,"\t",$Event[$i],"\n";
}

close(OUTPUT);

open OUTPUT,">$Target\_class.txt";

for (my $i=0;$i<@EventClass ;$i++) {
	print OUTPUT $EventClass[$i],"\n";
}

close(OUTPUT);


##########################################################################################################################################


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
