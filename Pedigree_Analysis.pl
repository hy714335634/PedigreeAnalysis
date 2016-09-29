use 5.010;				#版本声明;
use strict;				#严格语法;
use warnings;			#错误提示;
use File::Basename;		#路径操作;
use Cwd;				#路径操作;
#use utf8;				#中文;
#use open ":encoding(gbk)";	#中文;
use vars qw(%hash $Duration $os $Script_Start $cwd $path $Par $header);
$os = $^O ;			#系统判别;
$Script_Start=time;	#运行时间;
$header='Sample_ID	Index	Chromosome Position	Chr	Reference Nucleotide	Genotype	Coverage	A(#F,#R)	C(#F,#R)	G(#F,#R)	T(#F,#R)	Ins(#F,#R)	Del(#F,#R)	Low coverage flag	Strand bias flag	Multi-allele flag	Dense SNPs flag	comments	same_mutation_number	report_num	effect	Clinical significance	Disease	Inheritance	Main Clinical Manifestations	Gene	Mutation_Call:Relative to CDS	RNA_accession	mutation_type	Amino Acid Change	SNP db_xref	Mutant Allele Frequency	CDS	1000Gp1_ASN_AF	1000Gp1_AF	ExAC_AF	ExAC_EAS_AF	exome_AF	exome_EAS_AF	Minimum base	Minimum frequency	Mutation base	Mutation frequency	HGMD_ID	Gent_type	type	Variant_class	Disease	SIFT_score	PolyPhen2_HDIV_pred	PolyPhen2_HVAR_pred	Mutation_Taster_pred	PhyloP	CLINVAR_CLNSIG	CLINVAR_TRAIT';
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
my $temp_1;
my $Genotype_Select;
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



sub Help{						#帮助文档初始化;
	print "\n----------------------------Pedigree_Analysis.pl---------------------------\n"
		 ."\tFunction:Pedigree Analysis For WES Family\n"
		 ."\tUseage:'cd' the dir and run 'Perl Pedigree_Analysis.pl <-Parameter>'\n"
		 ."\tParameter:\n"
		 ."\t\t<-i>\tSample ID like '*_WES_Mutation_Report_output.txt'\n"
		 ."\t\t<-p>\tPaternal Sample ID\n"
		 ."\t\t<-m>\tMaternal Sample ID\n"
		 ."\t\tdefault:\ The Script Will Scan The Dir and Ensure ID\n"
		 ."\tOutput:\n\t\tCompound_Heterozygous.txt\n"
		 ."\t\tDenovo_Mutation.txt\n"
		 ."\t\tPaternal_Genetic.txt\n"
		 ."\t\tMaternal_Genetic.txt\n"
		 ."\n------------------------------------END-----------------------------------\n"
		 ."\t\tAuthor:ZQ\tDATA:2016-8-22\tVersion:1.0\n";
}

sub OpenDIR{			#根据系统获得当前目录地址;
	if ($0 =~ m{^/}) {			
	  $cwd = dirname($0);
	} else {
	  $cwd = dirname(getcwd()."/$0");
	}
	$path=$cwd;				#目录路径;
	if($os eq "MSWin32"){
		$$path=~s/\\/\\\\/g;			#根据系统选择斜杠和反斜杠;
	}
	return $path;
}

sub ExitFun{				#脚本结束初始化;
	$Duration = (time - $Script_Start)/60;						#统计时间;
	printf("\nThis Perl Script Has Been Running For %0.2f Minute\n",$Duration);
}
sub Analysis{				#分析流程;
	my($Sample,$Paternal,$Maternal)=@_;
	print "-----------------------------------------------------------------\n";
	print "Sample ID:$Sample\nPaternal ID:$Paternal\nMaternal ID:$Maternal\n";
	print "Check.";
	$|=1;
	open SAM,'<',"$Sample" or die "STOP\n>>Sample File $!\n";
	print ".........";
	foreach $line(<SAM>){					#读取SAM构建HASH;
		chomp($line);
		if($line=~/^[A-Z|\#|\s+]/){
			next;
		}
		@temp=split /\t/,$line;
		if(!($temp[5])){next;}
		push @CH,$line;				#存储样本数据分析复合杂合;
		if($temp[5]=~/;/){
			@temp_Genotype_1=split /;/,$temp[5];
			foreach $temp_gene(@temp_Genotype_1){
				if(($temp_gene=~/del/)||($temp_gene=~/ins/)){
					push @temp_Genotype,$temp_gene;
					next;
				}else{
					@temp_Genotype_2=split //,$temp_gene;
					push @temp_Genotype,@temp_Genotype_2;
				}
			}
		}else{
			if(($temp[5]=~/del/)||($temp[5]=~/ins/)){
				push @temp_Genotype,$temp[5];
			}else{
				@temp_Genotype=split //,$temp[5];
			}
		}
		foreach $temp(@temp_Genotype){
			chomp($temp);
			$temp=$temp[3].",".$temp[2].",".$temp;
			$Sample{$temp}=$line;
		}
	}
	close SAM;
	open PAT,'<',"$Paternal" or die "STOP\n>>Paternal File $!\n";
	print "................";
	foreach $line(<PAT>){					#读取PAT构建HASH;
		chomp($line);
		if($line=~/^[A-Z|\#|\s+]/){
			next;
		}
		@temp=split /\t/,$line;
		if(!($temp[5])){next;}
		if($temp[5]=~/;/){
			@temp_Genotype_1=split /;/,$temp[5];
			foreach $temp_gene(@temp_Genotype_1){
				if(($temp_gene=~/del/)||($temp_gene=~/ins/)){
					push @temp_Genotype,$temp_gene;
					next;
				}else{
					@temp_Genotype_2=split //,$temp_gene;
					push @temp_Genotype,@temp_Genotype_2;
				}
			}
		}else{
			if(($temp[5]=~/del/)||($temp[5]=~/ins/)){
				push @temp_Genotype,$temp[5];
			}else{
				@temp_Genotype=split //,$temp[5];
			}
		}
		foreach $temp(@temp_Genotype){
			chomp($temp);
			$temp=$temp[3].",".$temp[2].",".$temp;
			$Paternal{$temp}=$line;
		}
	}
	close PAT;
	open MAT,'<',"$Maternal" or die "STOP\n>>Maternal File $!\n";
	print "..........";
	foreach $line(<MAT>){					#读取MAT构建HASH;
		chomp($line);
		if($line=~/^[A-Z|\#|\s+]/){
			next;
		}
		@temp=split /\t/,$line;
		if(!($temp[5])){next;}
		if($temp[5]=~/;/){
			@temp_Genotype_1=split /;/,$temp[5];
			foreach $temp_gene(@temp_Genotype_1){
				if(($temp_gene=~/del/)||($temp_gene=~/ins/)){
					push @temp_Genotype,$temp_gene;
					next;
				}else{
					@temp_Genotype_2=split //,$temp_gene;
					push @temp_Genotype,@temp_Genotype_2;
				}
			}
		}else{
			if(($temp[5]=~/del/)||($temp[5]=~/ins/)){
				push @temp_Genotype,$temp[5];
			}else{
				@temp_Genotype=split //,$temp[5];
			}
		}
		foreach $temp(@temp_Genotype){
			chomp($temp);
			$temp=$temp[3].",".$temp[2].",".$temp;
			$Maternal{$temp}=$line;
		}
	}
	close PAT;
	#@Sample=<SAM>;@Paternal=<PAT>;@Maternal=<MAT>;
	print ".............";
	print "..........[OK]\n";
	foreach $line(sort keys %Sample){		#分类器，获得数据集;
		if(exists $Paternal{$line}){
			chomp($Sample{$line});
			$PG{$line}=$Sample{$line}."\n";
			#@temp=split /\t/,$Sample{$line};
			# if(($temp[5]=~/del/)||($temp[5]=~/ins/)){
				# #print "$temp[31]\n";
					# if($temp[31] ge 90){
						# $Gene{$temp[25]}=1;
					# }
			# }else{
				# @temp=split //,$temp[5];
				# next if(!($temp[1]));
				# if($temp[0]=$temp[1]){
					# @temp=split /\t/,$Sample{$line};
					# $Gene{$temp[25]}=1;
				# }else{
					# my @temp_Genotype=split /,/,$line;
					# my $var_1=$temp_Genotype[0].",".$temp_Genotype[1].",".$temp[0];
					# my $var_2=$temp_Genotype[0].",".$temp_Genotype[1].",".$temp[1];
					# if(((exists $Paternal{$var_1})&&(exists $Maternal{$var_2}))||((exists $Paternal{$var_2})&&(exists $Maternal{$var_1}))){
						# @temp=split /\t/,$Sample{$line};
						# $Gene{$temp[25]}=1;
					# }
				# }
			# }
		}
		if(exists $Maternal{$line}){
			chomp($Sample{$line});
			$MG{$line}=$Sample{$line}."\n";
			#@temp=split /\t/,$Sample{$line};
			# if(($temp[5]=~/del/)||($temp[5]=~/ins/)){
				# #print "$temp[31]\n";
					# if($temp[31] ge 90){
						# $Gene{$temp[25]}=1;
					# }
			# }else{
				# @temp=split //,$temp[5];
				# if($temp[0]=$temp[1]){
					# @temp=split /\t/,$Sample{$line};
					# $Gene{$temp[25]}=1;
				# }else{
					# my @temp_Genotype=split /,/,$line;
					# my $var_1=$temp_Genotype[0].",".$temp_Genotype[1].",".$temp[0];
					# my $var_2=$temp_Genotype[0].",".$temp_Genotype[1].",".$temp[1];
					# if(((exists $Paternal{$var_1})&&(exists $Maternal{$var_2}))||((exists $Paternal{$var_2})&&(exists $Maternal{$var_1}))){
						# @temp=split /\t/,$Sample{$line};
						# $Gene{$temp[25]}=1;
					# }
				# }
			# }
		}
		if(!(exists $Maternal{$line})&&!(exists $Paternal{$line})){			#新生突变;
			chomp($Sample{$line});
			$DM{$line}=$Sample{$line}."\n";
		}
	}
	#通过基因列表构建复合杂和数据集;
	# foreach $line(sort keys %Sample){
		# chomp($line);
		# chomp($Sample{$line});
		# @temp=split /\t/,$Sample{$line};
		# if($temp[25]){
			# $CH{$line}=$Sample{$line}."\n";
		# }
	# }
	#反转HASH去重复;
	%PG=reverse %PG;
	@PG=keys %PG;
	%PG=();
	# %CH=reverse %CH;
	# @CH=keys %CH;
	# %CH=();
	%DM=reverse %DM;
	@DM=keys %DM;
	%DM=();
	%MG=reverse %MG;
	@MG=keys %MG;
	%MG=();
	#重新构建索引;
	foreach $line(@PG){
		chomp($line);
		@temp=split /\t/,$line;
		$PG{$temp[1]}=$line;
	}
	foreach $line(@MG){
		chomp($line);
		@temp=split /\t/,$line;
		$MG{$temp[1]}=$line;
	}
	# foreach $line(@CH){
		# chomp($line);
		# @temp=split /\t/,$line;
		# $CH{$temp[1]}=$line;
	# }
	foreach $line(@DM){
		chomp($line);
		@temp=split /\t/,$line;
		$DM{$temp[1]}=$line;
	}
	#数据表头、正文输出;
	# open CH,'>',"Compound_Heterozygous\.txt";
	# print CH"$header\n";
	open DM,'>',"Denovo_Mutation\.txt";
	print DM"$header\n";
	open PG,'>',"Paternal_Genetic\.txt";
	print PG"$header\n";
	open MG,'>',"Maternal_Genetic\.txt";
	print MG"$header\n";
	foreach(sort{$a <=> $b}keys %PG){
		chomp($_);
		print PG"$PG{$_}\n";
	}
	print ">>Paternal_Genetic..............................................[OK]\n";
	foreach(sort{$a <=> $b}keys %MG){
		chomp($_);
		print MG"$MG{$_}\n";
	}
	print ">>Maternal_Genetic..............................................[OK]\n";
	foreach(sort{$a <=> $b}keys %DM){
		chomp($_);
		print DM"$DM{$_}\n";
	}
	print ">>Denovo_Mutation...............................................[OK]\n";
	# foreach(sort{$a <=> $b}keys %CH){
		# chomp($_);
		# push @CH_res,$CH{$_};
	# }
	#print ">>Compound_Heterozygous....................";
	close DM;close PG;close MG;
	$flag=1;
	my $i=0;
	my $last_gene="Understand";
	open CH,'>',"Compound_Heterozygous\.txt";
	print CH"$header\n";
	print ">>Compound_Heterozygous....................";
	foreach $line(@CH){				#复合杂合分析;
		chomp($line);
		if($line=~/^[A-Z|\#|\s+]/){next;}
		@temp=split /\t/,$line;
		if((!($temp[5]))||($temp[5]=~/-/)||($temp[5]=~/;/)||(!$temp[25])||($temp[31]=~/;/)||($temp[31] eq "")){					#去除不符合要求的点或不需要分析的点;
			next;
		}else{							   #待分析点的分析;
			if($last_gene eq "Understand"){			#确定第一个基因名;
				$last_gene=$temp[25];
				$flag=1;
				@temp_Genotype_1=();
			}
			if(($temp[25] eq $last_gene)&&($flag==1)){		#分析同一个基因的不同位点;
				if((($temp[5]=~/del/)||($temp[5]=~/ins/))&&($temp[31]>1)){	#考虑插入缺失的纯合分析,杂合不考虑;
					my $Freq=$temp[31]+0;
					if($Freq>=90){			#比率高于90%考虑，否则该基因不属于复合杂合;
						$temp=$temp[3].",".$temp[2].",".$temp[5];
						$temp_1=$temp[3].",".$temp[2].",".$temp[5];
					}else{
						$flag=0;
					}
					$Genotype_Select="Heterozygous";
				}else{					#常规点分析;
					@temp_Genotype_2=split //,$temp[5];
					$temp=$temp[3].",".$temp[2].",".$temp_Genotype_2[0];
					$temp_1=$temp[3].",".$temp[2].",".$temp_Genotype_2[1];
				}
				if($temp eq $temp_1){	#纯合点;
					if((exists $Maternal{$temp})||(exists $Paternal{$temp})){
						push @temp_Genotype_1,$line;
					}else{
						$flag=0;
					}
					$Genotype_Select="Homozygous";
				}else{					#杂合点要求分别在父母一方;
					if(((exists $Paternal{$temp})&&(exists $Maternal{$temp_1}))||((exists $Maternal{$temp})&&(exists $Paternal{$temp_1}))){
						push @temp_Genotype_1,$line;
					}else{
						$flag=0;
					}
					$Genotype_Select="Heterozygous";
				}
			}
			if($temp[25] ne $last_gene){			#下一个基因分析;
				if($flag==1){						#上一个基因属于复合杂合，输出到文件，初始化结果数组;
					if(@temp_Genotype_1==1){			#分析只有一个点的基因，除去杂合点基因;
						if($Genotype_Select eq "Homozygous"){
							$temp=shift @temp_Genotype_1;
							print CH"$temp\n";
						}
					}else{
						foreach $temp(@temp_Genotype_1){
							chomp($temp);
							print CH"$temp\n";
						}
					}
				}
				@temp_Genotype_1=();
				$last_gene=$temp[25];
				$flag=1;
				if((($temp[5]=~/del/)||($temp[5]=~/ins/))&&($temp[31]>1)){ #考虑插入缺失的纯合分析,杂合不考虑;
					my $Freq=$temp[31]+0;
					if($Freq >= 90){   #比率高于90%考虑，否则该基因不属于复合杂合;
						$temp=$temp[3].",".$temp[2].",".$temp[5];
						$temp_1=$temp[3].",".$temp[2].",".$temp[5];
					}else{
						$flag=0;
					}
					$Genotype_Select="Homozygous";
				}else{		#常规点分析;
					@temp_Genotype_2=split //,$temp[5];
					$temp=$temp[3].",".$temp[2].",".$temp_Genotype_2[0];
					$temp_1=$temp[3].",".$temp[2].",".$temp_Genotype_2[1];
				}
				if($temp eq $temp_1){		#纯合点;
					if((exists $Maternal{$temp})||(exists $Paternal{$temp})){
						push @temp_Genotype_1,$line;
					}else{
						$flag=0;
					}
					$Genotype_Select="Homozygous";
				}else{		#杂合点要求分别在父母一方;
					if(((exists $Paternal{$temp})&&(exists $Maternal{$temp_1}))||((exists $Maternal{$temp})&&(exists $Paternal{$temp_1}))){
						push @temp_Genotype_1,$line;
					}else{
						$flag=0;
					}
					$Genotype_Select="Heterozygous";
				}
			}
		}
	}
	print ".....................[OK]\n";
	close CH;
}
