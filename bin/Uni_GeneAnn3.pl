###     此程序进一步对重叠基因基于单个基因的分歧度进行过滤
###     useage :  perl Uni_GeneAnn2.pl All_uni_DB_adj_div.txt All_uni_DB_adj.pep 0.2

use strict;
use warnings;
use Bio::Perl;

my $DivFile = $ARGV[0];
my $PepFile = $ARGV[1];

my @b=split(/[\/\.]/,$PepFile); my $Target=$b[$#b-1];
my $Gene_Modeller=$ENV{'Gene_Modeller'};
my $Threshold = 0.2;
if ( $ARGV[2] ) {
	$Threshold = $ARGV[2];
}

my %DList = ();
open INPUT,"$DivFile";
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_); 
	$DList{$a[0]} = 1-$a[1];
}
close(INPUT);

my %SList = ();  my %NList = ();

open OUTPUT,">GeneLoc.txt";

my $seqdata = "$PepFile";
my $seqio   = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	if ( $DList{$seq->display_id} ) {
		my @a = split(/\s/,$seq->description);
		my @b = split(/\_/,$a[1]);
		my $Chr = "";  my $Beg = 0;  my $End = 0;
		if ( $a[0] =~ /(.*)\_(\d+)\_(\d+)$/) {
			$Chr = $1;
			$Beg = $2+$b[0]-1;
			$End = $2+$b[1]-1;
		}
		print OUTPUT $seq->display_id,"\t",$Chr,"\t",$Beg,"\t",$End,"\t",$a[2],"\t",$DList{$seq->display_id},"\n";
		$SList{$seq->display_id} = $seq->seq();
		$NList{$seq->display_id} = $Chr."_".$Beg."_".$End." ".$a[2];
	}
}	

close(OUTPUT);

`perl $Gene_Modeller/Tab_Rank.pl GeneLoc.txt GeneLoc_Sort.txt 2 3 4`;
`perl $Gene_Modeller/LocRegionUni.pl GeneLoc_Sort.txt $Threshold`;

my %KList =();
open INPUT,"GeneLoc_Sort_uni.txt";
open OUTPUT,">$Target\_Final.pep";

while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_); 
	$KList{$a[0]}+=1;
	print OUTPUT ">",$a[0]," ",$NList{$a[0]},"\n";
	print OUTPUT $SList{$a[0]},"\n";

}
close(INPUT);
close(OUTPUT);

open INPUT,"$DivFile";
open OUTPUT,">$Target\_Final.txt";

while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ( $KList{$a[0]} ) {
		print OUTPUT $_,"\n";
	}
}

close(INPUT);
close(OUTPUT);

########################################################################################################################################

my %AList = ();
open INPUT,"GeneLoc_Sort_class.txt";
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_);
	push @{$AList{$a[0]}},$a[1];
}
close(INPUT);

foreach my $key (keys %AList) {
	my @a = @{$AList{$key}};
	if ( $#a > 1 ) {
		for (my $i=0;$i<@a ;$i++) {
			if ( !$KList{$a[$i]} and $DList{$a[$i]} > 0.85 ) {
				print $a[$i],"\t",$DList{$a[$i]},"\n";
			}
		}
	}
}


unlink("GeneLoc.txt");
unlink("GeneLoc_Sort.txt");
unlink("GeneLoc_Sort_uni.txt");
