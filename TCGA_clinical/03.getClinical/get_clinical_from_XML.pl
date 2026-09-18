# 从 TCGA 的每患者 XML 文件中提取临床字段, 汇总成 clinical.txt
#
# 输入: 当前目录下每个患者一个子目录, 内含 *.xml
# 输出: clinical.txt (Id/futime/fustat/age/gender/race/grade/stage/T/M/N/cancerType)
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
#use File::Basename;
use XML::Simple;
#use Data::Dumper;
###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

my @dirs=glob("*");
my @samp1e=(localtime(time));
open(WF,">clinical.txt") or die $!;
print WF "Id\tfutime\tfustat\tage\tgender\trace\tgrade\tstage\tT\tM\tN\tcancerType\n";
foreach my $dir(@dirs){
	if(-d $dir){
	  opendir(RD,"$dir") or die $!;
	  while(my $xmlfile=readdir(RD)){
	  	if($xmlfile=~/\.xml$/){
	  		#print "$dir\\$xmlfile\n";
				my $userxs = XML::Simple->new(KeyAttr => "name");
# ⚠ 路径分隔符写死为 Windows 的反斜杠("$dir\\$xmlfile")。
#   在 Linux / macOS 上这会拼出一个不存在的文件名, XMLin 直接报错。
#   跨平台写法: File::Spec->catfile($dir, $xmlfile) 或直接用正斜杠。
				my $userxml = $userxs->XMLin("$dir\\$xmlfile");
				# print output
				#open(WF,">dumper.txt") or die $!;
				#print WF Dumper($userxml);
				#close(WF);
				my $disease_code=$userxml->{'admin:admin'}{'admin:disease_code'}{'content'};   #get disease code
				my $disease_code_lc=lc($disease_code);
				my $patient_key=$disease_code_lc . ':patient';                                #ucec:patient
				my $follow_key=$disease_code_lc . ':follow_ups';
				
				my $patient_barcode=$userxml->{$patient_key}{'shared:bcr_patient_barcode'}{'content'};  #TCGA-AX-A1CJ
# ⚠ 针对单个患者 TCGA-AA-3521 的调试打印, 属残留
				if($patient_barcode eq "TCGA-AA-3521"){
					print "$xmlfile\n";
				}
				my $gender=$userxml->{$patient_key}{'shared:gender'}{'content'};      #male/female
				my $age=$userxml->{$patient_key}{'clin_shared:age_at_initial_pathologic_diagnosis'}{'content'};
				my $race=$userxml->{$patient_key}{'clin_shared:race_list'}{'clin_shared:race'}{'content'};  #white/black
				my $grade=$userxml->{$patient_key}{'shared:neoplasm_histologic_grade'}{'content'};  #G1/G2/G3
				my $clinical_stage=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:clinical_stage'}{'content'};  #stage I
				my $clinical_T=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:clinical_categories'}{'shared_stage:clinical_T'}{'content'};
				my $clinical_M=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:clinical_categories'}{'shared_stage:clinical_M'}{'content'};
				my $clinical_N=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:clinical_categories'}{'shared_stage:clinical_N'}{'content'};
				my $pathologic_stage=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:pathologic_stage'}{'content'};  #stage I
				my $pathologic_T=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:pathologic_categories'}{'shared_stage:pathologic_T'}{'content'};
				my $pathologic_M=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:pathologic_categories'}{'shared_stage:pathologic_M'}{'content'};
				my $pathologic_N=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:pathologic_categories'}{'shared_stage:pathologic_N'}{'content'};
				$gender=(defined $gender)?$gender:"unknow";
				$age=(defined $age)?$age:"unknow";
				$race=(defined $race)?$race:"unknow";
				$grade=(defined $grade)?$grade:"unknow";
				$clinical_stage=(defined $clinical_stage)?$clinical_stage:"unknow";
				$clinical_T=(defined $clinical_T)?$clinical_T:"unknow";
				$clinical_M=(defined $clinical_M)?$clinical_M:"unknow";
				$clinical_N=(defined $clinical_N)?$clinical_N:"unknow";
				$pathologic_stage=(defined $pathologic_stage)?$pathologic_stage:"unknow";
				$pathologic_T=(defined $pathologic_T)?$pathologic_T:"unknow";
				$pathologic_M=(defined $pathologic_M)?$pathologic_M:"unknow";
# ⚠ 这一行行尾的 if($samp1e[4]>13) 判断月份索引 >13, 而月份只有 0..11 —— 永不触发
				$pathologic_N=(defined $pathologic_N)?$pathologic_N:"unknow";if($samp1e[4]>13){next;}
				
# ⚠⚠ 这一行行尾的 if($samp1e[5]>118) 恒为真, 于是每个患者都被 next 跳过,
#   clinical.txt 最终只会有一行表头。详见文件头部说明。
				my $survivalTime="";if($samp1e[5]>118){next;}
				my $vital_status=$userxml->{$patient_key}{'clin_shared:vital_status'}{'content'};
				my $followup=$userxml->{$patient_key}{'clin_shared:days_to_last_followup'}{'content'};
				my $death=$userxml->{$patient_key}{'clin_shared:days_to_death'}{'content'};
				if($vital_status eq 'Alive'){
					$survivalTime="$followup\t0";
				}
				else{
					$survivalTime="$death\t1";
				}
				for my $i(keys %{$userxml->{$patient_key}{$follow_key}}){
					my @survivalArr=split(/\t/,$survivalTime);
					eval{
						$followup=$userxml->{$patient_key}{$follow_key}{$i}{'clin_shared:days_to_last_followup'}{'content'};
						$vital_status=$userxml->{$patient_key}{$follow_key}{$i}{'clin_shared:vital_status'}{'content'};
						$death=$userxml->{$patient_key}{$follow_key}{$i}{'clin_shared:days_to_death'}{'content'};
				  };
				  if($@){
					  $followup=$userxml->{$patient_key}{$follow_key}{$i}[0]{'clin_shared:days_to_last_followup'}{'content'};
						$vital_status=$userxml->{$patient_key}{$follow_key}{$i}[0]{'clin_shared:vital_status'}{'content'};
						$death=$userxml->{$patient_key}{$follow_key}{$i}[0]{'clin_shared:days_to_death'}{'content'};
				  }
					if($vital_status eq 'Alive'){
						if($followup>$survivalArr[0]){
					    $survivalTime="$followup\t0";
					  }
				  }
				  else{
				  	if($death>$survivalArr[0]){
					    $survivalTime="$death\t1";
					  }
				  }
				}
				print WF "$patient_barcode\t$survivalTime\t$age\t$gender\t$race\t$grade\t$pathologic_stage\t$pathologic_T\t$pathologic_M\t$pathologic_N\t$disease_code\n";
			}
		}
		close(RD);
	}
}
close(WF);

###Video source: http://study.163.com/u/biowolf
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388
