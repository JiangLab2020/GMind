use strict;
use warnings;

my $input     =    $ARGV[0];             ####  记录排序前的文件
my $output    =    $ARGV[1];             ####  记录排序后的文件

`dos2unix $input`;

my $SortFlag1 =    1;                    ####  第一排序列, 默认对第一列进行排序
if ( $ARGV[2] ) {
	$SortFlag1 =   $ARGV[2];
}
my $SortFlag2 =    "";                   ####  第二排序列
if ( $ARGV[3] ) {
	$SortFlag2 =   $ARGV[3];
}

my $SortFlag3 =    "";                   ####  第三排序列
if ( $ARGV[4] ) {
	$SortFlag3 =   $ARGV[4];
}

my @Mark = ();                           ####  记录每一列是字符串(1)还是数字(2)

my $reg1 = qr/^-?\d+(\.\d+)?$/;          ####  匹配数字
my $reg2 = qr/^-?0(\d+)?$/;

open INPUT,"$input";
my @line=();
while (<INPUT>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ( @a ) {
		for (my $i=0;$i<@a ;$i++) {
			if ( !$Mark[$i] or $Mark[$i] == 2 ) {
				if ( $a[$i] eq 0 or ($a[$i] =~ $reg1 && $a[$i] !~ $reg2) ) {
					$Mark[$i] = 2;
				}
				else {
					$Mark[$i] = 1;
				}
			}
		}	
		push @line,[@a];
	}
}
close(INPUT);

open OUTPUT,">$output";
if ( @line ) {
	my @linesort=();
	if ( $SortFlag3 ) {                                  ###     必须保证第一个是字符,第三个是数字
		$SortFlag1-=1; $SortFlag2-=1; $SortFlag3-=1;
		print $SortFlag1,"\t",$Mark[$SortFlag1],"\n";
		print $SortFlag2,"\t",$Mark[$SortFlag2],"\n";
		if ( $Mark[$SortFlag2] == 2 ) {
			@linesort=sort { $a->[$SortFlag1] cmp $b->[$SortFlag1] or $a->[$SortFlag2] <=> $b->[$SortFlag2] or $a->[$SortFlag3] <=> $b->[$SortFlag3]  } @line;
		}
		else {
			@linesort=sort { $a->[$SortFlag1] cmp $b->[$SortFlag1] or $a->[$SortFlag2] cmp $b->[$SortFlag2] or $a->[$SortFlag3] <=> $b->[$SortFlag3]  } @line;
		}
	}
	elsif ( $SortFlag2 ) {
		my $Col1 = abs($SortFlag1)-1;
		my $Col2 = abs($SortFlag2)-1;
		if ( $SortFlag1 > 0 and $SortFlag2 > 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $a->[$Col1] cmp $b->[$Col1] or $a->[$Col2] cmp $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 > 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $a->[$Col1] cmp $b->[$Col1] or $a->[$Col2] <=> $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 > 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $a->[$Col1] <=> $b->[$Col1] or $a->[$Col2] cmp $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 > 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $a->[$Col1] <=> $b->[$Col1] or $a->[$Col2] <=> $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 < 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $b->[$Col1] cmp $a->[$Col1] or $b->[$Col2] cmp $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 < 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $b->[$Col1] cmp $a->[$Col1] or $b->[$Col2] <=> $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 < 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $b->[$Col1] <=> $a->[$Col1] or $b->[$Col2] cmp $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 < 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $b->[$Col1] <=> $a->[$Col1] or $b->[$Col2] <=> $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 < 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $a->[$Col1] cmp $b->[$Col1] or $b->[$Col2] cmp $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 < 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $a->[$Col1] cmp $b->[$Col1] or $b->[$Col2] <=> $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 < 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $a->[$Col1] <=> $b->[$Col1] or $b->[$Col2] cmp $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 > 0 and $SortFlag2 < 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $a->[$Col1] <=> $b->[$Col1] or $b->[$Col2] <=> $a->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 > 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $b->[$Col1] cmp $a->[$Col1] or $a->[$Col2] cmp $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 > 0 and $Mark[$Col1] == 1 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $b->[$Col1] cmp $a->[$Col1] or $a->[$Col2] <=> $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 > 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 1 ) {
			@linesort=sort { $b->[$Col1] <=> $a->[$Col1] or $a->[$Col2] cmp $b->[$Col2] } @line;
		}
		elsif ( $SortFlag1 < 0 and $SortFlag2 > 0 and $Mark[$Col1] == 2 and $Mark[$Col2] == 2 ) {
			@linesort=sort { $b->[$Col1] <=> $a->[$Col1] or $a->[$Col2] <=> $b->[$Col2] } @line;
		}
	}
	else {
		my $Col1 = abs($SortFlag1)-1;
		if ( $SortFlag1 > 0 and $Mark[$Col1] == 1 ) {
			@linesort=sort { $a->[$Col1] cmp $b->[$Col1] } @line;
		}
		elsif ( $SortFlag1 > 0 and $Mark[$Col1] == 2 ) {
			@linesort=sort { $a->[$Col1] <=> $b->[$Col1] } @line;
		}
		elsif ( $SortFlag1 < 0 and $Mark[$Col1] == 1 ) {
			@linesort=sort { $b->[$Col1] cmp $a->[$Col1] } @line;
		}
		elsif ( $SortFlag1 < 0 and $Mark[$Col1] == 2 ) {
			@linesort=sort { $b->[$Col1] <=> $a->[$Col1] } @line;
		}
	}
	foreach ( @linesort )
	{
		print OUTPUT join("\t",@$_),"\n";
	}
}
close(OUTPUT);
