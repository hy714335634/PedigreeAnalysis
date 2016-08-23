use 5.010;				#版本声明;
use strict;				#严格语法;
use warnings;			#错误提示;
use File::Basename;		#路径操作;
use Cwd;				#路径操作;
use utf8;				#中文;
use open ":encoding(gbk)";	#中文;
my %hash;				#哈希结构;
my $Duration;			#时间统计; 
my $os = $^O ;			#系统判别;
my $Script_Start=time;	#运行时间;
my $cwd;				#路径变量;
my $path;				#路径变量;
my $Par;				#参数下标;
my $header='Sample_ID	Index	Chromosome Position	Chr	Reference Nucleotide	Genotype	Coverage	A(#F,#R)	C(#F,#R)	G(#F,#R)	T(#F,#R)	Ins(#F,#R)	Del(#F,#R)	Low coverage flag	Strand bias flag	Multi-allele flag	Dense SNPs flag	comments	same_mutation_number	report_num	effect	Clinical significance	Disease	Inheritance	Main Clinical Manifestations	Gene	Mutation_Call:Relative to CDS	RNA_accession	mutation_type	Amino Acid Change	SNP db_xref	Mutant Allele Frequency	CDS	1000Gp1_ASN_AF	1000Gp1_AF	ExAC_AF	ExAC_EAS_AF	exome_AF	exome_EAS_AF	Minimum base	Minimum frequency	Mutation base	Mutation frequency	HGMD_ID	Gent_type	type	Variant_class	Disease	SIFT_score	PolyPhen2_HDIV_pred	PolyPhen2_HVAR_pred	Mutation_Taster_pred	PhyloP	CLINVAR_CLNSIG	CLINVAR_TRAIT';

#-------------------------脚本变量-----------------------#
my $Sample_name;		#样本名;
my $Paternal_name;		#父系样本名;
my $Maternal_name;		#母系样本名;
my $flag=0;				#是否符合运行要求;
my $file;				#临时文件;
my @temp_Genotype;
my @temp_Genotype_1;
my @temp_Genotype_2;
#--------------------------------------------------------#
my %Sample;				#主样本;
my %Paternal;			#父样本;
my %Maternal;			#母样本;
my %CH;					#存储复合杂合;
my @CH;					#存储复合杂合;
my %DM;					#存储新生突变;
my @DM;					#存储新生突变;
my %PG;					#存储父源遗传;
my @PG;					#存储父源遗传;
my %MG;					#存储母源遗传;
my @MG;					#存储母源遗传;
my %Gene;				#存储复合杂合相关基因;
my $temp;				#临时变量;
my @temp;				#临时存储;
my $line;				#临时变量;
my $temp_gene;			#临时变量;
#------------------------参数配置------------------------#
my @Parameter=@ARGV;
if(@Parameter){					#获得参数;
	foreach $Par(0..$#Parameter){
		given($Parameter[$Par]){				#读取参数;
				when($_ eq '-i') {
				$Sample_name=$Parameter[$Par+1];
				$flag+=1;
				break;
			}
			when($_ eq '-p') {
				$Paternal_name=$Parameter[$Par+1];
				$flag+=1;
				break;
			}
			when($_ eq '-m') {
				$Maternal_name=$Parameter[$Par+1];
				$flag+=1;
				break;
			}
			when(($_ eq '-help')||($_ eq 'help\/\?')) {
				&Help;
				$flag=-1;
				break;
			}
			default {
				if($Parameter[$Par]=~/^-/){		#默认;
					&Help;
					$flag=-1;
					last;
				}
			}
		}
	}
}
if($flag==3){
	&Analysis($Sample_name,$Paternal_name,$Maternal_name);
	&ExitFun;
	$flag=-1;
}elsif($flag==0){
	$path=&OpenDIR;
	opendir(DIR,"$path");			#打开当前目录;
	while($file=readdir(DIR)){
		if(!($file=~/^\.+$/)){
			given($file){				#自动筛选文件;
				when($_=~/^\d+_WES_Mutation/){
					$Sample_name=$file;
					$flag+=1;
					break;
				}
				when($_=~/^\d+\-0201+_WES_Mutation/){
					$Paternal_name=$file;
					$flag+=1;
					break;
				}
				when($_=~/^\d+\-0301+_WES_Mutation/){
					$Maternal_name=$file;
					$flag+=1;
					break;
				}
			}
		}
	}
	if($flag==3){				#判断是否满足执行条件;
		&Analysis($Sample_name,$Paternal_name,$Maternal_name);
		&ExitFun;
	}else{
		print "\n>>Can't Find Correct File as Sample,Paternal,Maternal!\n";
	}
}elsif($flag>0){
	print "\n>>Please Input 'Perl Pedigree_Analysis.pl <-i> Parameter <-p> Parameter <-m> Parameter.\n>>Or You can run 'Perl Pedigree_Analysis.pl -help' For Help !\n";
}
