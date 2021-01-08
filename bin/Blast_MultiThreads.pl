#!/usr/bin/perl -w

####   此程序用于多线程基因HMM

use strict;
use warnings;

my $database =     "$ARGV[0]";             ####  需要查询物种的Fasta文件
my $query    =     "$ARGV[1]";             ####  需要查询物种的Fasta文件
my $outputTem=     "$ARGV[2]";             ####  此文件用于记录每次一个物种对一个物种Blast的结果,
my $AlignParameter="$ARGV[3]";             ####  表示是进行哪种比对,如果参数值是1,则表示蛋白质比对,
                                           ####                    ;如果参数值是2,则表示核酸比对.
								           ####                    ;如果参数值是3,则表示核酸到蛋白质的比对.
										   ####                    ;如果参数值是4,则表示蛋白质到核酸的比对.

my $EValue=1e-6;                          ####  局部比对的阈值

my @b=split(/[\/\.]/,$query); my $species=$b[$#b-1]."_Tem";

my $ReNum    = 1;                          ####  同一个查询所能找到的靶基因的最大值
if ( $ARGV[4] ) {
	$ReNum=$ARGV[4];
}

my $Group = 1;                             ####  CPU数量.
if ( $ARGV[5] ) {
	$Group=$ARGV[5];
}

my $MakeDB = 1;                            ####  是否建库,默认建库.
if ( $ARGV[6] ) {
	$MakeDB=$ARGV[6];
}


print "####  1 Counting the number of sequences\n";
########################################################################################################
my $Len=0;  my @AllGeneName = ();
open INPUT,$query;
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
	$Group=$Len;
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
open INPUT,$query;
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
my $formatdb="makeblastdb";
my $blastall=""; my $dbtype=""; my $suff="";


if ( $AlignParameter == 1) {
	$blastall="blastp";
	$dbtype="prot";
	$suff="p";
}
if ( $AlignParameter == 2) {
	$blastall="blastn";
	$dbtype="nucl";
	$suff="n";
}
if ( $AlignParameter == 3) {
	$blastall="blastx";
	$dbtype="prot";
	$suff="p";
}
if ( $AlignParameter == 4) {
	$blastall="tblastn";
	$dbtype="nucl";
	$suff="n";
}

if ( $database eq $query and $ReNum == 1 ) {
	$ReNum=50;
}

#if ( -e "/public/home/chengj/Programs/blast/bin" ) {
#	$formatdb="/public/home/chengj/Programs/blast/bin/$formatdb";
#	$blastall="/public/home/chengj/Programs/blast/bin/$blastall";
#}

if ( $MakeDB == 1 ) {
	`$formatdb -in $database -dbtype $dbtype -parse_seqids -max_file_sz '3GB'`; 
}

if ( -e "$database."."$suff"."to" and -e "$database."."$suff"."db" ) {
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
			
			`$blastall -db $database -query $species\_$j\.fas -out $species\_$j\.txt -max_target_seqs $ReNum -evalue $EValue -outfmt 6`;		

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

	open OUTPUT,">$outputTem";

	for (my $j=0;$j<$Group ;$j++) {	
		open INPUT,"$species\_$j\.txt";
		while (<INPUT>) {
			chomp($_);
			print OUTPUT $_,"\n";
		}
		close(INPUT);
		unlink("$species\_$j\.txt");
		unlink("$species\_$j\.fas");
	}

	close(OUTPUT);


	########################################################################################################
	if ( $MakeDB == 1 ) {
		unlink("$database."."$suff"."hr");
		unlink("$database."."$suff"."in");
		unlink("$database."."$suff"."si");
		unlink("$database."."$suff"."sd");
		unlink("$database."."$suff"."sq");
		unlink("$database."."$suff"."og");
		unlink("$database."."$suff"."nd");
		unlink("$database."."$suff"."ni");

		unlink("$database."."$suff"."db");
		unlink("$database."."$suff"."os");
		unlink("$database."."$suff"."ot");
		unlink("$database."."$suff"."tf");
		unlink("$database."."$suff"."to");
	}
}
else {
	print "Error in makeblastdb!!!\n";
}
