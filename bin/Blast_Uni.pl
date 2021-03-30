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

my $GMind=$ENV{'GMind'};
my $GMind_lib=$ENV{'GMind_lib'};
`perl $GMind/Tab_Rank.pl GeneLocPS.txt GeneLocPS_Sort.txt 1 3 4`;
`perl $GMind/Tab_Rank.pl GeneLocNS.txt GeneLocNS_Sort.txt 1 3 4`;

`perl $GMind/LocRegionUni.pl GeneLocPS_Sort.txt $Threshold`;
`perl $GMind/LocRegionUni.pl GeneLocNS_Sort.txt $Threshold`;

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

unlink("GeneLocPS.txt");
unlink("GeneLocPS_Sort.txt");
unlink("GeneLocPS_Sort_uni.txt");
unlink("GeneLocPS_Sort_class.txt");

unlink("GeneLocNS.txt");
unlink("GeneLocNS_Sort.txt");
unlink("GeneLocNS_Sort_uni.txt");
unlink("GeneLocNS_Sort_class.txt");
