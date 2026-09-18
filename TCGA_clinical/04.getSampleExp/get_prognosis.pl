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
			if(($samples[3]=~/^1/)&&($normalFlag==1)){
				push(@indexs,$i);
				push(@normalSamples,$arr[$i]);
			}
			else{
				my $sampleName="$samples[0]-$samples[1]-$samples[2]";
				if(exists $hash1{$sampleName}){
					push(@indexs1,$i);
					push(@samples1,$arr[$i]);if($samp1e[5]>118){next;}
					delete($hash1{$sampleName});
				}
				if(exists $hash2{$sampleName}){
					push(@indexs2,$i);
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