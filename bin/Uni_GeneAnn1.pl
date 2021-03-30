use strict;
use warnings;
use Bio::Perl;

my @AnnFiles = ();
if ( $ARGV[0] ) {
	push(@AnnFiles,$ARGV[0]);
}
if ( $ARGV[1] ) {
	push(@AnnFiles,$ARGV[1]);
}
if ( $ARGV[2] ) {
	push(@AnnFiles,$ARGV[2]);
}
if ( $ARGV[3] ) {
	push(@AnnFiles,$ARGV[3]);
}
my $GMind=$ENV{'GMind'};

my $Query_DB = "Query_DB.pep";

my $BlastNum = 20;

my $CPU_Num = 20;

my %SList = ();  my %NList = ();

my $seqio   = new Bio::SeqIO(-file => $Query_DB, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	$SList{$seq->display_id} = $seq->seq();
}	

my %KList =();

open OUTPUT,">GeneLoc.txt";

for (my $i=0;$i<@AnnFiles ;$i++) {
	my $seqdata = $AnnFiles[$i];
	my $seqio   = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
	while ( my $seq = $seqio->next_seq ) {
		my $seqlength=length($seq->seq());
		my $MySeq = $seq->seq();   $MySeq =~ s/\*$//;
		
		my $Conut = 0; $Conut=$MySeq=~s/\*/\*/g; $Conut++;  $KList{$seq->display_id} = $Conut;
		my $Score = 500/$Conut-abs($seqlength-500);
		my @a = split(/\s/,$seq->description);
		my @b = split(/\_/,$a[1]);
		print OUTPUT $seq->display_id,"\t",$a[0],"\t",$b[0],"\t",$b[1],"\t",$a[2],"\t",$Score,"\n";
		$SList{$seq->display_id} = $seq->seq();
		$NList{$seq->display_id} = $seq->description;
	}	
}

close(OUTPUT);

`perl $GMind/Tab_Rank.pl GeneLoc.txt GeneLoc_Sort.txt 2 3 4`;
`perl $GMind/LocRegionUni.pl GeneLoc_Sort.txt`;


open INPUT,"GeneLoc_Sort_uni.txt";
open OUTPUT,">GeneLoc_Sort_uni.pep";

while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_); 
	print OUTPUT ">",$a[0],"\n";
	print OUTPUT $SList{$a[0]},"\n";
}

close(INPUT);
close(OUTPUT);

if ( !(-e "GeneLoc_Sort_uni_DB.txt") ) {
	`perl $GMind/Blast_MultiThreads.pl $Query_DB GeneLoc_Sort_uni.pep GeneLoc_Sort_uni_DB.txt 1 $BlastNum $CPU_Num`;
}

my %BList = ();   my %TList = ();
open INPUT,"GeneLoc_Sort_uni_DB.txt";
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_);
	if ( !$TList{$a[0].$a[1]} ) {
		$TList{$a[0].$a[1]}+=1;
		push @{$BList{$a[0]}},$a[1];
	}
}
close(INPUT);


my %AList = ();
open INPUT,"GeneLoc_Sort_class.txt";
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_);
	push @{$AList{$a[0]}},$a[1];
}
close(INPUT);

my %DIvst=();
open OUTPUT1,">GeneLoc_div.txt";
open OUTPUT,">GeneLoc_Sort_uni_DB.pep";

foreach my $key (keys %AList) {
	my @a = @{$AList{$key}};
	print $key,"\n";

	my %PList = (); my %MList = ();

	my $AlignNum = 0;

	open OUT,">Tem_Tem.pep";
	for (my $i=0;$i<@a ;$i++) {
		$PList{$a[$i]}+=1;
		print OUT ">",$a[$i],"\n",$SList{$a[$i]},"\n";
		$AlignNum++;
		if ( $a[$i]=~ /\_(.*)\_/ ) {
			push @{$MList{$1}},$a[$i];
		}
		if ( $BList{$a[$i]} ) {
			my @b = @{$BList{$a[$i]}};
			for (my $j=0;$j<@b ;$j++) {
				print OUT ">",$b[$j],"\n",$SList{$b[$j]},"\n";
				$AlignNum++;
			}
		}
	}
	close(OUT);

	if ( $AlignNum > 1 ) {	
		`clustalo -i Tem_Tem.pep -o Tem_Tem.fas -v --outfmt fa --threads 10 --force --output-order tree-order`;
		`perl $GMind/Fasta_DelLargeDiv.pl Tem_Tem.fas Tem_Tem.txt 0`;

		my $MinName = "";  my $MinDivg = 1;
		open IN,"Tem_Tem.txt";
		while (<IN>) {
			chomp($_);
			my @a = split(/\t/,$_);
			$DIvst{$a[0]}=$a[1];
			if ( $PList{$a[0]} and $a[1] < $MinDivg and $KList{$a[0]} < 3 ) {
				$MinName = $a[0];  $MinDivg = $a[1];
			}
		}
		close(IN);

		if ( $MinName ) {
			if ( $MinName =~ /\_(.*)\_/ ) {
				my @a = @{$MList{$1}};
				for (my $j=0;$j<@a ;$j++) {
					print OUTPUT ">",$a[$j]," ",$NList{$a[$j]},"\n",$SList{$a[$j]},"\n";
					print OUTPUT1 $a[$j],"\t",$DIvst{$a[$j]},"\n";
				}
			}
		}
	}
}
close(OUTPUT);
close OUTPUT1;

unlink("GeneLoc.txt");
unlink("GeneLoc_Sort.txt");
unlink("GeneLoc_Sort_uni_DB.txt");
unlink("GeneLoc_Sort_uni.txt");
unlink("GeneLoc_Sort_uni.pep");
