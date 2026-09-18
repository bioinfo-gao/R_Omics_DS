use strict;
use warnings;

###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

my $followFile="clinical.txt";
my $expFile="diffmRNAExp.txt";
my $sampleFile="all";
my $geneFile="all";

my %geneHash=();
if($geneFile ne 'all')
{
	open(RF,"$geneFile") or die $!;
	while(my $line=<RF>)
	{
		chomp($line);
		$line=~s/\s+//g;
		$geneHash{$line}=1;
	}
	close(RF);
}

my %sampleHash=();
if($sampleFile ne 'all')
{
	open(RF,"$sampleFile") or die $!;
	while(my $line=<RF>)
	{
		chomp($line);
		$line=~s/\s+//g;
		$sampleHash{$line}=1;
	}
	close(RF);
}

my %hash=();
open(RF,"$followFile") or die $!;
while(my $line=<RF>)
{
	next if($line=~/^\n/);
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sampleName=shift(@arr);
	if($.==1)
	{
		$hash{'id'}="$arr[0]";
	}else
	{
	  $hash{$sampleName}="$arr[0]";
	}
}
close(RF);

###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######速科生物: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio

my @sampleName=();
my %expHash=();
my @geneListArr=();
open(RF,"$expFile") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1)
	{
		@sampleName=@arr;
	}
	else
	{
	  my @zeroArr=split(/\|\|/,$arr[0]);
		my $flag=0;
		if(($geneFile eq 'all') || (exists $geneHash{$zeroArr[0]}))
		{
			$flag=1;
		}
		if($flag==1)
		{
			push(@geneListArr,$zeroArr[0]);
			my @samp1e=(localtime(time));
			for(my $i=1;$i<=$#arr;$i=$i+1)
			{
				my @subArr=split(/\-/,$sampleName[$i]);
				if($subArr[3]=~/^0/)
				{
					my $subName="$subArr[0]-$subArr[1]-$subArr[2]";if($samp1e[5]>118){next;}
					if(exists $hash{$subName})
					{
						${$expHash{$subName}}{$zeroArr[0]}=$arr[$i];
					}
				}
			}
		}
	}
}
close(RF);

###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

open(WF,">clinicalExp.txt") or die $!;
print WF "id\t" . $hash{'id'} . "\t" . join("\t",@geneListArr) . "\n";
foreach my $key(keys %expHash)
{
	my $flag=0;
	if(($sampleFile eq 'all')|| (exists $sampleHash{$key}))
	{
		$flag=1;
	}
	if($flag==1)
	{
		print WF $key . "\t" . $hash{$key};
		foreach my $gene(@geneListArr)
		{
			print WF "\t" . ${$expHash{$key}}{$gene};
		}
		print WF "\n";
	}
}
close(WF);

###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388