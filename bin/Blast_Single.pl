#!/usr/bin/perl -w

# 此程序用于做简单Blast。
use strict;
use warnings;
use Bio::Perl;
###############################################################################################     参数设定

my $database=      "$ARGV[0]";             ####  需要查询物种的Fasta文件
my $query =        "$ARGV[1]";             ####  需要查询物种的Fasta文件
my $outputTem=     "$ARGV[2]";             ####  此文件用于记录每次一个物种对一个物种Blast的结果,
my $AlignParameter="$ARGV[3]";             ####  表示是进行哪种比对,如果参数值是1,则表示蛋白质比对,
                                           ####                    ;如果参数值是2,则表示核酸比对.
								           ####                    ;如果参数值是3,则表示核酸到蛋白质的比对.
										   ####                    ;如果参数值是4,则表示蛋白质到核酸的比对.
my $EValue=1e-6;                          ####  局部比对的阈值

my $ReNum    = 1;                          ####  同一个查询所能找到的靶基因的最大值
if ( $ARGV[4] ) {
	$ReNum=$ARGV[4];
} 

my $Identity = 0;                          ####  仅取出Identity大于该阈值的比对.
if ( $ARGV[5] ) {
	$Identity=$ARGV[5];
} 

#print OUTPUT1 "Query id,Subject id,% identity,alignment length,mismatches,gap openings,q.start,q.end,s.start,s.end,e-value,bit score,P(A),Species\n";

#################################################################################################### Blast搜索						

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

################################################################################################## Blast搜索

#if ( -e "/public/home/chengj/Programs/blast/bin" ) {
#	$formatdb="/public/home/chengj/Programs/blast/bin/$formatdb";
#	$blastall="/public/home/chengj/Programs/blast/bin/$blastall";
#}


`$formatdb -in $database -dbtype $dbtype -parse_seqids -max_file_sz '3GB'`;

if ( $ReNum > 0 ) {
	`$blastall -db $database -query $query -out $outputTem -max_target_seqs $ReNum -evalue $EValue -outfmt 6`;
}
elsif ( $ReNum == -1 ) {
	`$blastall -db $database -query $query -out $outputTem -num_threads 4 -evalue $EValue -outfmt 6`;     ###  最多输出500条序列
}
else {
	`$blastall -db $database -query $query -out $outputTem -num_threads 4 -max_target_seqs 10000 -evalue $EValue -outfmt 6`;
}


#`$blastall -db $database -query $query -out $outputTem -evalue $EValue`;


############################################################################################################

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

if ( $Identity > 0 ) {
	my %LList   =();
	my %PList   =();
	if ( -e "$outputTem" ) {
		my $seqdata = $query;
		my $seqio   = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
		while ( my $seq = $seqio->next_seq ) {
			my $MySeq = $seq->seq();  $MySeq =~ s/\*$//;
			$LList{$seq->display_id}=length($MySeq);
		}
	}
	if ( -e "$outputTem" ) {
		my $seqdata = $database;
		my $seqio   = new Bio::SeqIO(-file => $seqdata, -format => 'fasta');
		while ( my $seq = $seqio->next_seq ) {
			my $MySeq = $seq->seq();  $MySeq =~ s/\*$//;
			$PList{$seq->display_id}=length($MySeq);
		}
	}	

	`cp $outputTem $outputTem\.txt`;

	open INPUT1,"$outputTem\.txt"; 
	open OUTPUT1,">$outputTem";
	while (<INPUT1>) {
		chomp($_);
		my @a=split(/\t/,$_);
		if ( !$PList{$a[1]} ) {
			print $a[1],"\n";
			$a[1]=lc($a[1]);
		}
		if ( $a[2] > $Identity ) {
			my $StrTem = join("\t",@a);
			print OUTPUT1 $StrTem,"\t",$LList{$a[0]},"\t",Float($a[3]/$LList{$a[0]}),"\t",$PList{$a[1]},"\t",Float($a[3]/$PList{$a[1]}),"\n";
		}	
	}
	close(INPUT1);
	close(OUTPUT1);

	unlink("$outputTem\.txt");
}




############    保留小数点后两位
sub Float{
	my $str=shift;
	if ($str =~ /\./) {
		$str = sprintf "%0.2f",$str;
	}
	return $str;
}
