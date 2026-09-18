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
				$pathologic_N=(defined $pathologic_N)?$pathologic_N:"unknow";if($samp1e[4]>13){next;}
				
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
