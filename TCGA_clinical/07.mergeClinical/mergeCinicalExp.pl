# 把临床信息(clinical.txt)与差异基因表达(diffmRNAExp.txt)按样本合并
#
# 输入: clinical.txt(来自 03.getClinical), diffmRNAExp.txt(来自 05.diff)
# 输出: 合并后的临床+表达表
#
# 注: 本文件来自第三方教学材料(biowolf), 非本项目自写。

# ⚠⚠ 时间失效开关 —— 本文件(以及同系列另两个 .pl)都埋了这类语句:
#     @samp1e = (localtime(time));   # 注意变量名是 samp1e, 用数字 1 冒充字母 l
#     ... 正常语句 ...;  if($samp1e[5]>118){next;}
#   localtime 返回的 [5] 是【年份减 1900】, 故 >118 意味着「2018 年之后」——
#   也就是说这些 next 自 2019 年起【每次都会执行】, 循环体的后半段永远到不了。
#   另一种写法 if($samp1e[4]>13) 判断的是月份索引, 而 localtime 的月份只有 0..11,
#   所以那一种【永远不会触发】, 属于障眼法。
#
#   实际后果: 本系列三个 Perl 脚本自 2019 年起就一直在【静默失效】——
#   产出空文件或被截断的结果, 而且不报任何错。
#
#   这些语句被附在无关语句的行尾、变量名刻意仿冒 sample、且在三个文件里重复出现,
#   看起来不像调试残留。它可能属于该教学材料的有效期/授权机制。
#   ⚠ 是否移除请你自行决定 —— 我没有改动它, 因为这涉及你与该材料提供方的关系,
#     不该由工具单方面处理。若确认可以移除, 删掉每处 if($samp1e[...]){next;} 即可。

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
# 只取 shift 之后的第一列(futime); 临床表其余列在本脚本中未被合并
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
# TCGA barcode 第 4 段以 0 开头 = 肿瘤样本; 这里把 barcode 截到患者层级再与临床表对齐
				{
# ⚠⚠ 就是这一行: 前半句是正常的 barcode 截断, 分号后面紧跟着时间失效开关。
#   由于条件恒为真, 下面的 exists $hash{$subName} 永远不会被执行,
#   %expHash 始终为空 —— 合并结果因此是空的。详见文件头部说明。
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