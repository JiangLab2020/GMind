use strict;
use warnings;

my $InputFA =  $ARGV[0];                     ####  需要库的Fasta文件

my $Group = 1;                               ####  需要将查询的Fasta文件分成多少组 
if ( $ARGV[1] ) {
	$Group  =  $ARGV[1];
}

my $Gene_Modeller_lib=$ENV{'Gene_Modeller_lib'};
my $PfamPath = "$Gene_Modeller_lib\/Pfam-A.hmm";

my $PF_ID   = 0;                             ####  默认是注释所有家族, 或者指定家族PF_ID.
if ( $ARGV[2] ) {
	$PF_ID  = $ARGV[2];
	$PfamPath = "$Gene_Modeller_lib\/HMM/$PF_ID\.hmm";
	if ( !(-e "$Gene_Modeller_lib\/HMM/$PF_ID\.hmm.h3f") ) {
		if ( -e "$Gene_Modeller_lib\/HMM/$PF_ID\.hmm.h3i" ) {
			unlink("$Gene_Modeller_lib\/HMM/$PF_ID\.hmm.h3i");
		}
		`hmmpress $PfamPath`;
	}	
}

my $All_Flag = 0;                            ####  是否输出所有结果.
if ( $ARGV[3] ) {
	$All_Flag = $ARGV[3];
}

my @b=split(/[\/\.]/,$InputFA); my $species=$b[$#b-1];

my $PValue   =  1e-4;


print "####  1 Counting the number of sequences\n";
########################################################################################################
my $Len=0;  my @AllGeneName = ();
open INPUT,$InputFA;
while (<INPUT>) {
	if ($_ =~ /^\>(.*)/) {
		my @a=split(/\s+/,$1);
		push @AllGeneName,$a[0];
		$Len++;
	}
}
close(INPUT);
print "       ","Len = ","$Len","\n";
if ( $Group > $Len ) {
	$Group=1;
}
########################################################################################################

print "####  2 Spliting the sequences to $Group groups\n";
########################################################################################################
my @GroLim=();  my @Cut=();
my $GroLen=int $Len/$Group;
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

print "####  3 MultiThreads Hmmer\n";
########################################################################################################
my @childs=();
my $count=0;

for (my $j=0;$j<$Group ;$j++) {
	$count++;                          #   print "fork process $count\n";
	my $pid = fork();			
	if ($pid) {                        #   parent
		sleep 1;
		push(@childs, $pid);
	}
	elsif ($pid == 0) {
		print "       ","fork process $count\n";
		
		`hmmscan --cpu 10 --domtblout $species\_$j\.txt -E $PValue $PfamPath $species\_$j\.fas`;

		exit(0);
	}
	else {
		die "couldn't fork: $!\n";
	}
}
foreach (@childs) {
	waitpid($_,0);
}
########################################################################################################

print "####  4 Fetching protein sequences from domtblout file\n";
########################################################################################################

my %CList = ();  my %DList = ();  my %EList = ();

open OUTPUT,">$species\_Hmm_All.txt";

for (my $j=0;$j<$Group ;$j++) {	
	open INPUT,"$species\_$j\.txt";
	while (<INPUT>) {
		chomp($_);
		if ( $_ =~ /^\#/ ) {
			next;
		}
		print OUTPUT $_,"\n";
		my @a=split(/\s+/,$_);

		my @b=(); $b[0]=$a[3]; $b[1]=$a[1]; $b[2]=$a[11]; $b[3]=$a[17]; $b[4]=$a[18];

		push @{$EList{$a[3]}},[@b];

		if ( !$DList{$a[1]."\t".$a[3]} ) {
			$DList{$a[1]."\t".$a[3]}+=1;
			$CList{$a[3]}.="\t".$a[1];
		}
	}
	close(INPUT);
	unlink("$species\_$j\.txt");
	unlink("$species\_$j\.fas");
}

close(OUTPUT);

open OUTPUT,">$species\_Hmm.txt";

foreach  (@AllGeneName) {
	if ( $CList{$_} ) {
		print OUTPUT $_,$CList{$_},"\n";
	}
	else {
		print OUTPUT $_,"\n";
	}
}

close(OUTPUT);

if ( $All_Flag == 0 ) {
	unlink("$species\_Hmm_All.txt");
}

if ( $All_Flag == 1 ) {
	my %FList = ();

	open OUTPUT,">$species\_Hmm_Uni.txt";

	foreach my $key (keys %EList) {
		my @line = @{$EList{$key}};
		my @linesort=sort { $a->[2] <=> $b->[2] } @line;
		my @Uni = ();
		foreach ( @linesort ) {
			my @a = @$_; my $Flag = 0;
			if ( @Uni ) {
				for (my $i=0;$i<@Uni ;$i++) {
					my @b = @{$Uni[$i]};
					my $OverlapLen = RegionOverlap("$b[3]\t$b[4]","$a[3]\t$a[4]");
					if ( $OverlapLen/($a[4]-$a[3]+1) > 0.8 ) {
						$Flag = 1;  last;
					}
				}
			}
			if ( $Flag == 0 ) {
				push @Uni,[@a];
			}
		}

		@linesort=sort { $a->[3] <=> $b->[3] } @Uni;
		foreach ( @linesort ) {
			my @a = @$_;
			push @{$FList{$a[0]}},$a[1];
		}
	}

	foreach my $key (sort keys %FList) {
		print OUTPUT $key,"\t",join("\t",@{$FList{$key}}),"\n";
	}

	close(OUTPUT);
}

########################################################################################################



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

