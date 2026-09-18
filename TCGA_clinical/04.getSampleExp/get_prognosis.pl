# 按 sample1.txt / sample2.txt 指定的两组样本, 从 input.txt 抽取表达矩阵
#
# 输入: input.txt(表达矩阵), sample1.txt, sample2.txt(各一行一个样本名)
# 输出: sampleExp.txt (供 05.diff/edgeR.R 使用)
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
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

my %hash1=();
my %hash2=();
my $normalFlag=0;
my $sampleFile1="sample1.txt";
my $sampleFile2="sample2.txt";

open(RF,"$sampleFile1") or die $!;
while(my $line=<RF>){
	chomp($line);
	$hash1{$line}=1;
}
close(RF);
my @samp1e=(localtime(time));
open(RF,"$sampleFile2") or die $!;
while(my $line=<RF>){
	chomp($line);
	$hash2{$line}=1;
}
close(RF);

my @indexs=();
my @indexs1=();
my @indexs2=();
open(RF,"input.txt") or die $!;
open(WF,">sampleExp.txt") or die $!;
my @normalSamples=();
my @samples1=();
my @samples2=();
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		for(my $i=1;$i<=$#arr;$i++){
			my @samples=split(/\-/,$arr[$i]);
# normalFlag 默认为 0, 故第 4 段以 1 开头的正常样本这一支【从不进入】——
# 想同时取正常样本需把上面的 $normalFlag 改成 1
			if(($samples[3]=~/^1/)&&($normalFlag==1)){
				push(@indexs,$i);
				push(@normalSamples,$arr[$i]);
			}
			else{
				my $sampleName="$samples[0]-$samples[1]-$samples[2]";
				if(exists $hash1{$sampleName}){
					push(@indexs1,$i);
# ⚠⚠ 行尾的 if($samp1e[5]>118) 恒为真, 会在表头解析途中跳出, 使样本选取被截断
					push(@samples1,$arr[$i]);if($samp1e[5]>118){next;}
					delete($hash1{$sampleName});
				}
				if(exists $hash2{$sampleName}){
					push(@indexs2,$i);
# ⚠ 行尾的 if($samp1e[4]>13) 判断月份 >13, 永不触发
					push(@samples2,$arr[$i]);if($samp1e[4]>13){next;}
					delete($hash2{$sampleName});
				}
			}
		}
		if($normalFlag==1){
			print WF "ID\t" . join("\t",@normalSamples) . "\t" . join("\t",@samples1) . "\t" . join("\t",@samples2) . "\n";
		}
		else{
		  print WF "ID\t" . join("\t",@samples1) . "\t" . join("\t",@samples2) . "\n";
	  }
	}
	else{
		my @sampleNor=();
		my @sampleData1=();
		my @sampleData2=();
		if($normalFlag==1){
		  foreach my $col(@indexs){
			  push(@sampleNor,$arr[$col]);
		  }
	  }
		foreach my $col(@indexs1){
			push(@sampleData1,$arr[$col]);
		}
		foreach my $col(@indexs2){
			push(@sampleData2,$arr[$col]);
		}
		if($normalFlag==1){
			print WF "$arr[0]\t" . join("\t",@sampleNor) . "\t" . join("\t",@sampleData1) . "\t" . join("\t",@sampleData2) . "\n";
		}
		else{
		  print WF "$arr[0]\t" . join("\t",@sampleData1) . "\t" . join("\t",@sampleData2) . "\n";
	  }
	}
}
close(WF);
close(RF);

#print "normal: " . ($#normalSamples+1) . "\n";
print "sample1: " . ($#samples1+1) . "\n";
print "sample2: " . ($#samples2+1) . "\n";

###Video source: http://study.163.com/u/biowolf
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388