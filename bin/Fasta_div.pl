use strict;
use warnings;
use Bio::Perl;

my $InputFA =  $ARGV[0];                     ####  输入需计算分歧度的Fasta文件
my $QueryDB =  $ARGV[1];                     ####  需要库的Fasta文件
my $Group   =  $ARGV[2];                     ####  需要将查询的Fasta文件分成多少组  

my @b=split(/[\/\.]/,$InputFA); my $species=$b[$#b-1];
my $Gene_Modeller=$ENV{'Gene_Modeller'};
my $Gene_Modeller_lib=$ENV{'Gene_Modeller_lib'};

my %SList=();								####  记录库中序列
my $seqdata1 = $QueryDB;
my $seqio1 = new Bio::SeqIO(-file => $seqdata1, -format => 'fasta');
while ( my $seq = $seqio1->next_seq ) {   
	$SList{$seq->display_id}=$seq->seq();
}
print "####  1 Counting the number of sequences\n";
########################################################################################################
my $Len=0;
my %AList = ();             ###  记录每条序列
my %BList = ();				###  记录物种名
my $seqdata = $InputFA;
my $seqio = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {   
	$AList{$seq->display_id}=$seq->seq();
	$BList{$seq->display_id}=$seq->description;
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

print "####  3 MultiThreads Blast\n";
########################################################################################################
my $BlastNum   = 20;

my @childs=();
my $count=0;

open OUTPUT3,">$species\_div.txt";
for (my $j=0;$j<$Group ;$j++) {
	$count++;                          #   print "fork process $count\n";
	`cp $QueryDB $QueryDB\_$j`;
	my $pid = fork();			
	if ($pid) {                        #   parent
		sleep 1;
		push(@childs, $pid);
	}
	elsif ($pid == 0) {
		print "       ","fork process $count\n";
		
		`perl $Gene_Modeller/Blast_Single.pl $QueryDB\_$j $species\_$j\.fas $species\_$j\.txt 1 $BlastNum`;

		my %TList=();
		my %CList=();
		open INPUT,"$species\_$j\.txt";
		while (<INPUT>) {
			chomp($_);
			my @a=split(/\t/,$_);
			$CList{$a[0]}=1;
			push @{$TList{$a[0]}},$a[1];
		}
		close INPUT;
		foreach my $key (sort keys %CList) {

			open OUT,">Tem_$j\.pep";
			print OUT ">",$key,"\n",$AList{$key},"\n";
			my @b=@{$TList{$key}};
			for (my $k=0;$k<@b ;$k++) {
				if ( !$CList{$b[$k]} ) {
					print OUT ">",$b[$k],"\n",$SList{$b[$k]},"\n";
				}
			}
			close OUT;

			`clustalo -i Tem_$j\.pep -o Tem_$j\.fas -v --outfmt fa --threads 2 --force --output-order tree-order`;
			`perl $Gene_Modeller/Fasta_DelLargeDiv.pl Tem_$j\.fas Tem_$j\.txt 0`;

			open IN,"Tem_$j\.txt";
			while (<IN>) {
				chomp($_);
				my @d=split(/\t/,$_);
				if ($CList{$d[0]}) {
					print OUTPUT3 $d[0],"\t",$d[1],"\n";
				}
			}
			close IN;
		}
	
		exit(0);
	}
	else {
		die "couldn't fork: $!\n";
	}
}
foreach (@childs) {
	waitpid($_,0);
}

close OUTPUT3;

for (my $j=0;$j<$Group ;$j++) {
	unlink("$species\_$j\.txt");
	unlink("$species\_$j\.fas");
	unlink("Tem_$j\.txt");
	unlink("Tem_$j\.fas");
	unlink("Tem_$j\.pep");
	unlink("$QueryDB\_$j");
}


########################################################################################################