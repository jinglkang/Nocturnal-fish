# Nocturnal fish analysis restart
## Core clock genes
```bash
# kangjingliang@KangdeMBP-2 日 12 21 2025 15:03:22 ~/Documents/2025/Nocturnal_fish/Specific_genes
perl extract_target_genes.pl Core_clock_genes.txt Core_clock_genes
cat *.fas > CoreCR.fasta
scp CoreCR.fasta jlkang@10.33.247.14:~/Nocturnal_fish/Orthologous/pep/CoreCR

# jlkang@hnu2024 Sun Jan 04 2026 16:54:30 ~/Nocturnal_fish/Orthologous/pep
mkdir CoreCR
cp ../Phototransduction/all_test.fas ./
# blast all sequences of the test species to Phototransduction.fasta
# diamond
# jlkang@hnu2024 Sun Jan 04 2026 17:04:49 ~/Nocturnal_fish/Orthologous/pep/CoreCR
diamond makedb --in CoreCR.fasta -d CoreCR
diamond blastp -q all_test.fas -e 1e-5 --sensitive -k 1 -d CoreCR --out CoreCR_blastp.result

# prepare the fasta file per gene for the input of orthofinder
# jlkang@hnu2024 Sun Jan 04 2026 17:07:21 ~/Nocturnal_fish/Orthologous/pep/CoreCR
cp ../Phototransduction/*.pl ./
perl prepare_seq_orthofinder.pl CoreCR_blastp.result CoreCR.fasta

# run orthofinder for each gene
for i in *_gene;do cd ${i};orthofinder -f ./;cd ../;done
# select the orthogroups with 80% nocturnal (23) and diurnal (8) fish species transcripts
perl select_qualified_orth.pl > Qualified_orth.txt
less Qualified_orth.txt|cut -f 1|sort -u|wc -l # 17: no duplicated genes

# obtain the orthologous_list_rep.txt (only keep the representative transcript per species)
perl create_orth_rep_list.pl CoreCR_blastp.result CoreCR.fasta
# prepare the input sequences for evolutionary analyse
perl create_orth_seq.pl CoreCR_blastp.result # all sequences are saved in "sequences/"

# create the "correlation.txt" in "*_gene/OrthoFinder/Results_Dec14/Orthogroups"
# Results_Jan04
perl create_correlation.pl
perl Recode.pl > Qualified_orth_recode.txt

# picked the required species for phylogenetic tree construction
# based on the previous orthologous genes protein alignment
# some phylogeny trees have been built in the "phototransduction"
# make sure which we should build
# jlkang@hnu2024 Sun Jan 04 2026 21:19:41 ~/Nocturnal_fish/Orthologous/pep/CoreCR
perl check_phy.pl > Qualified_orth_recode_check.txt
cp ../Phototransduction/conca_sigpep.fa.phy ./
# jlkang@hnu2024 Sun Jan 04 2026 21:31:14 ~/Nocturnal_fish/Orthologous/pep/CoreCR
nohup perl Build_phylogeny.pl > tree.reports 2>&1 &
# [1] 1293381

# cp tree to the target directory
# jlkang@hnu2024 Sun Jan 04 2026 23:44:09 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl cp_tree.pl

# paml: free-ratio
```

```check_phy.pl
#!/usr/bin/perl
use strict;
use warnings;

# Phototransduction/Qualified_orth_recode.txt
# CoreCR/Qualified_orth_recode.txt
my $pho="../Phototransduction/Qualified_orth_recode.txt";
my %hash;
open PHO, $pho or die "can not open $pho\n";
while (<PHO>) {
        chomp;
        my @a=split /\t/;ls
        my $nb=$a[2]+$a[3];
        $hash{$a[4]}="Orth_".$nb."spe_".$a[-1]."code";
}

my $ccr="Qualified_orth_recode.txt";
open CCR, $ccr or die "can not open $ccr\n";
while (<CCR>) {
        chomp;
        my @a=split;
        my $nb=$a[2]+$a[3];
        my $code="Orth_".$nb."spe_".$a[-1]."code";
        if ($hash{$a[4]}) {
                print "$_\t$code\t$hash{$a[4]}\n";
        } else {
                print "$_\t$code\n";
        }
}
```

```cp_tree.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# ~/Nocturnal_fish/Orthologous/pep/Phototransduction/cnga1a_gene/OrthoFinder/Results_Dec14/Orthogroups/OG0000000
my $code="Qualified_orth_recode.txt";
open CODE, $code or die "can not open $code\n";
while (<CODE>) {
	chomp;
	my @a=split /\t/;
	my $nb  =$a[-1];
	my $tree="spe_".$nb."code.tre";
	my $rep =$a[0]."_gene/OrthoFinder/Results_Dec14/Orthogroups/final_orth_input_paml.txt";
	open REP, $rep or die "can not open $rep\n";
	while (<REP>) {
		chomp;
		my $dir=dirname($rep);
		$dir=$dir."/".$_;
		my $cmd="cp Phylogeny\/$tree $dir\/spe.tre";
		#print "$cmd\n";
		system($cmd);
	}
}
```
## Use Orthofinder
```bash
# (base) shichuang001@login01 五 3月 20 2026 17:47:13 ~/jlkang/Nocturnal_fish/Orthologous/pep
mamba activate orthofinder_env

# h2076@h2076 Wed Apr 01 2026 19:58:50 ~/Nocturnal_fish/Opsins
conda config --add channels http://sxycloud1.top:11428/conda/ekernf01
conda config --add channels http://sxycloud1.top:11428/conda/radiomics
conda config --add channels http://sxycloud1.top:11428/conda/dacase
conda config --add channels http://sxycloud1.top:11428/conda/pyg
conda config --add channels http://sxycloud1.top:11428/conda/terhorst
conda config --add channels http://sxycloud1.top:11428/conda/r
conda config --add channels http://sxycloud1.top:11428/conda/msys2
conda config --add channels http://sxycloud1.top:11428/conda/pytorch
conda config --add channels http://sxycloud1.top:11428/conda/nvidia
conda config --add channels http://sxycloud1.top:11428/conda/conda-forge
conda config --add channels http://sxycloud1.top:11428/conda/bioconda
conda config --add channels http://sxycloud1.top:11428/conda/free
conda config --add channels http://sxycloud1.top:11428/conda/main
conda config --set show_channel_urls yes
# h2076@h2076 Wed Apr 01 2026 20:06:34 ~/Nocturnal_fish/Opsins
conda install orthofinder -c bioconda
# To activate this environment, use
conda activate orthofinder_env
# To deactivate an active environment, use
conda deactivate

# Select the orthogroups with at least 60% nocturnal (17) and diurnal (6) fish species which must include at least one zebrafish transcript
# h2076@h2076 Sat Apr 04 2026 18:09:32 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp2.pl > orthologous_list_rep_2.txt
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

# 60%
# 28 nocturnal species: at least 17 species (60%)
# 10 diurnal species: at least 6 species (60%)
# make sure the slected orthogroups including Zebrafish

my %durnal=(
    'Acura'=> 1,
    'Apoly'=> 1,
    'Daru'=> 1,
    'Pmol'=> 1,
    'Padel'=> 1,
    'Platyfish'=> 1,
    'Fugu'=> 1,
    'Medaka'=> 1,
    'Stickleback'=> 1,
    'Zebrafish'=> 1,
);


my @headers; my %hash1; my $j;
my $numb="Orthogroups.GeneCount.tsv";
open NUMB, $numb or die "can not open $numb\n";
while (<NUMB>) {
    chomp;
    s/\r//g;
    my @a=split /\t/;
    if (/^Orthogroup/) {
        @headers=@a;
        #print "$_\n";
    } else {
        my ($durn, $noct);
        for (my $i = 1; $i < @a-1; $i++) {
            if ($durnal{$headers[$i]} && $a[$i] >= 1) {
                $durn++;
            } elsif ($a[$i] >= 1) {
                $noct++;
            }
        }
        #print "$_\n" if ($durn && $noct) && ($durn >= 8 && $noct >= 23) && ! ($durn == 10 && $noct == 28);
        #print "$a[0]\t$durn\t$noct\n" if ($durn && $noct) && ($durn >= 8 && $noct >= 23) && ! ($durn == 10 && $noct == 28);
        my $info;
        if (($durn && $noct) && ($durn >= 6 && $noct >= 17) && ($a[-4] >= 1) && ! ($durn == 10 && $noct == 28)) {
            $hash1{$a[0]}++;
            # print "$_\n";
        }
    }
}

my @spes=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus 
    Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus 
    Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma 
    Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora 
    Zleptacanthus Zviridiventer Acura Apoly Daru Pmol Padel 
    Platyfish Fugu Medaka Stickleback Zebrafish);
my $head="Orth\t";
foreach my $spe (@spes) {
    $head.=$spe."\t";
}
$head=~s/\s+$//;
print "$head\n";

foreach my $orth (sort keys %hash1) {
    #print "$orth\n";
    my $fasta="../Orthogroup_Sequences/$orth.fa";
    my ($spe, $orthid, %orth);
    open FASTA, $fasta or die "can not open $fasta\n";
    while (<FASTA>) {
        chomp;
        if (/\>/) {
            s/\>//; $orthid=$_;
            ($spe)=$orthid=~/(.*)\_.*/;
        } else {
            if ($orth{$spe}) {
                my $oldlen=length($orth{$spe}->{'SEQ'});
                my $newlen=length($_);
                if ($oldlen < $newlen) {
                    $orth{$spe}={
                        'ID'  => $orthid,
                        'SEQ' => $_ 
                    };
                }
            } else {
                $orth{$spe}={
                    'ID'  => $orthid,
                    'SEQ' => $_ 
                };  
            }
        }
    }
    my $info="$orth\t";
    foreach my $sp (@spes) {
        my $ID;
        ($orth{$sp}->{'ID'})?($ID=$orth{$sp}->{'ID'}):($ID="None");
        $info.=$ID."\t";
    }
    $info=~s/\s+$//;
    print "$info\n";
}
```

```bash
# Annotation
# h2076@h2076 Sat Apr 04 2026 18:35:20 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp3.pl orthologous_list_rep_2.txt > orthologous_list_rep_2_ano.txt
perl temp3.pl orthologous_list_rep.txt > orthologous_list_rep_ano.txt
# revise temp3.pl for the input file
```

```temp3.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

my %hash;
my $rep=$ARGV[0];
open REP, $rep or die "can not open $rep\n";
while (<REP>) {
    chomp;
    next if /^Orth/;
    my @a=split /\t/;
    $hash{$a[-1]}=$a[0];
}

my $ano="all_swissprot_diamond_ano_final.txt";
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
    chomp;
    my @a=split /\t/;
    if ($hash{$a[0]}) {
        print "$hash{$a[0]}\t$a[-0]\t$a[1]\t$a[-1]\n";
    }
}
```

## Core clock genes
```bash
# h2076@h2076 Sun Apr 05 2026 08:47:53 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
vi Core_clock_genes.txt
# ROR
# NR1D
# NFIL3
# HLF
# TEF
# CLOCK
# BMAL
# CRY
# PER
# BHE40
# BHE41
# CIART
# DBP

# h2076@h2076 Sun Apr 05 2026 10:40:17 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp4.pl Core_clock_genes.txt orthologous_list_rep_ano.txt > CCGs_rep.txt
perl temp4.pl Core_clock_genes.txt orthologous_list_rep_2_ano.txt > CCGs_rep_2.txt
# Edit Orth_CCGs.txt to keep the required genes
```

```temp4.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

# Pick the core clock genes
# RORB, RORAA, NR1D1, NR1D2, NFIL3
# HLF, TEF, CLOCK, BMAL1, CRY, PER
# BHE40, BHE41, CIART, DBP

my %hash;
my $gene=$ARGV[0];
open GENE, $gene or die "can not open $gene\n";
while (<GENE>) {
    chomp;
    $hash{$_}++;
}

my $ano=$ARGV[1];
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
    chomp;
    my @a=split /\t/;
    my $orth;
    ($orth)=$a[2]=~/sp\|.*\|(.*?)_.*/;
    my $lenA=length($orth);
    foreach my $key (sort keys %hash) {
        my $lenB=length($key);
        if ($lenB > $lenA) {
            next;
        } else {
            my $sub;
            $sub=substr($orth, 0, $lenB);
            if ($sub eq $key) {
                print "$_\n";
                last;
            }
        }
    }
}
```

## Phototransduction genes
```bash
# phototransduction genes except opsins
# h2076@h2076 Sun Apr 05 2026 11:20:07 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
less Phototransduction_genes.txt|perl -alne 'my @a=split /\t/;(my $nm)=$a[1]=~/(.*)\_.*/;print "$nm\t$a[-1]"' > Phototransduction_genes_uniq.txt
# kangjingliang@KangdeMBP-2 日  4 05 2026 12:10:01 ~/Desktop
scp -r -P 20231 h2076@sxycloud1.top:~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/Phototransduction_genes_uniq.txt ./
scp -r -P 20231 Phototransduction_genes_uniq.txt h2076@sxycloud1.top:~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/
# h2076@h2076 Sun Apr 05 2026 12:20:01 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp5.pl Phototransduction_genes_uniq.txt orthologous_list_rep_ano.txt > Phototransduction_orth_genes_rep.txt
perl temp5.pl Phototransduction_genes_uniq.txt orthologous_list_rep_2_ano.txt > Phototransduction_orth_genes_rep2.txt
cat Phototransduction_orth_genes_rep.txt Phototransduction_orth_genes_rep2.txt|perl -alne 'my @a=split /\t/;($nm)=$a[2]=~/sp\|.*\|(.*)_.*/;$a[3]=~s/\s+\[.*\]$//;print "$nm\t$a[3]"'|sort -u
less Phototransduction_genes_uniq.txt|perl -alne 'print $F[0]' > Phototransduction_genes_uniq_name.txt
cat Phototransduction_genes_uniq_name.txt opsin_genes.txt > phototransduction_genes_total_name.txt
perl temp9.pl|sort -u > phototransduction_genes_total_name_target.txt

# opsins
# kangjingliang@KangdeMBP-2 日  4 05 2026 12:16:51 ~/Desktop
less uniprotkb_opsin_AND_reviewed_true_2026_04_05.tsv|perl -alne 'next if /^Entry/;my @a=split /\t/;my $nm;($nm)=$a[2]=~/(.*?)\_/;print $nm if $nm=~/^OP/'|sort -u > opsin_genes.txt
scp -r -P 20231 opsin_genes.txt h2076@sxycloud1.top:~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/
# h2076@h2076 Sun Apr 05 2026 12:20:37 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp4.pl opsin_genes.txt orthologous_list_rep_2_ano.txt > opsins_rep_2.txt # opsins in at least 17 nocturnal and 6 diurnal species (60% and include zebrafish)
perl temp4.pl opsin_genes.txt orthologous_list_rep_ano.txt > opsins_rep.txt # opsins in all 28 nocturnal and 10 diurnal species
```

```temp9.pl
#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $name="phototransduction_genes_total_name.txt";
open NAME, $name or die "can not open $name\n";
while (<NAME>) {
    chomp;
    $hash{$_}++;
}

my $anno="all_swissprot_diamond_ano_final.txt";
open ANNO, $anno or die "can not open $anno\n";
while (<ANNO>) {
    chomp;
    my @a=split /\t/;
    my $nm;
    ($nm)=$a[1]=~/sp\|.*\|(.*)\_.*/;
    $a[-1]=~s/\s+?\[.*\]$//;
    print "$nm\t$a[-1]\n" if $hash{$nm};
}
```

## paml analysis
```bash
# h2076@h2076 Sun Apr 05 2026 16:51:11 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
cat Orth_CCGs.txt opsins_rep_2.txt Phototransduction_orth_genes.txt > paml_rep_2_genes.txt

# check how many trees we need construct
# h2076@h2076 Sun Apr 05 2026 17:31:15 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp6.pl orthologous_list_rep_2.txt paml_rep_2_genes.txt > paml_rep_2_genes_treeInfo.txt
# 20 trees
# h2076@h2076 Sun Apr 05 2026 17:35:16 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
less conca_sigpep.fa # pick the corresponding species for tree construction
# h2076@h2076 Sun Apr 05 2026 21:15:50 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
cp conca_sigpep.fa ../

# prepare the sequences for phylogeny tree construction
# h2076@h2076 Sun Apr 05 2026 21:36:28 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp7.pl # all the fasta file were put into "trees/"

# transform the fasta file to phylip file
# h2076@h2076 Sun Apr 05 2026 22:05:33 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/trees
for i in *.fa; do name=$(basename ${i} .fa);fasta2phy.pl ${i} > ${name}.phy; done
# scp to another workstation for tree construction
# (base) shichuang001@login01 日 4月 05 2026 22:15:14 ~/jlkang/Nocturnal_fish
scp -r -P 20231 h2076@sxycloud1.top:~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/trees ./
# (base) shichuang001@login01 日 4月 05 2026 22:55:08 ~/jlkang/Nocturnal_fish/trees
sbatch tree1.job
# Submitted batch job 113798
sbatch tree2.job
# Submitted batch job 113799
sbatch tree3.job
# Submitted batch job 113803

# 好慢而且一晚上才跑了4个bootstraps，然后就到时间了
# 换到数信跑
# select the best model
# h2076@h2076 Mon Apr 06 2026 08:47:34 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/trees
nohup ProteinModelSelection.pl tree1.phy > model.outfile 2>&1 &
# [1] 824351
# Best Model : JTTF
```

```tree1.job
#!/bin/bash
#SBATCH -J tree1 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 192 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.tree1 # 标准输出文件
#SBATCH -e logs/%j.tree1 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s tree1.phy -o Zebrafish -n tree1 -T 192 # 执行命令
echo end on $(date) # 结束时间
```

```tree2_ModelSelection.job
#!/bin/bash
#SBATCH -J tree2Model # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.tree2Model # 标准输出文件
#SBATCH -e logs/%j.tree2Model # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree2.phy > tree2Model.outfile # 执行命令
echo end on $(date) # 结束时间
```

```bash
# (base) shichuang001@login01 一 4月 06 2026 10:45:02 ~/jlkang/Nocturnal_fish/trees
sbatch tree1.job
# Submitted batch job 113804
# (base) shichuang001@login01 一 4月 06 2026 10:50:50 ~/jlkang/Nocturnal_fish/trees
sbatch tree2_ModelSelection.job
# Submitted batch job 113805
# (base) shichuang001@login01 一 4月 06 2026 12:51:52 ~/jlkang/Nocturnal_fish/trees
cat tree2Model.outfile
# Best Model : JTTF
# 全部的Best Model都是JTTF
```

```tree2.job
#!/bin/bash
#SBATCH -J tree2 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 192 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.tree2 # 标准输出文件
#SBATCH -e logs/%j.tree2 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree2.phy -o Zebrafish -n tree2 -T 192 # 执行命令
echo end on $(date) # 结束时间
```

```tree_ModelSel_3-5.job
#!/bin/bash
#SBATCH -J treeM3-5 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeM3-5 # 标准输出文件
#SBATCH -e logs/%j.treeM3-5 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree3.phy > tree3Model.outfile # 执行命令
ProteinModelSelection.pl tree4.phy > tree4Model.outfile
ProteinModelSelection.pl tree5.phy > tree5Model.outfile
echo end on $(date) # 结束时间
```

```tree_ModelSel_6-8.job
#!/bin/bash
#SBATCH -J treeM6-8 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeM6-8 # 标准输出文件
#SBATCH -e logs/%j.treeM6-8 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree6.phy > tree6Model.outfile
ProteinModelSelection.pl tree7.phy > tree7Model.outfile
ProteinModelSelection.pl tree8.phy > tree8Model.outfile
echo end on $(date) # 结束时间
```

```tree_ModelSel_9-11.job
#!/bin/bash
#SBATCH -J treeM9-11 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeM9-11 # 标准输出文件
#SBATCH -e logs/%j.treeM9-11 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree9.phy > tree9Model.outfile
ProteinModelSelection.pl tree10.phy > tree10Model.outfile
ProteinModelSelection.pl tree11.phy > tree11Model.outfile
echo end on $(date) # 结束时间
```

```tree_ModelSel_12-14.job
#!/bin/bash
#SBATCH -J treeM12-14 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeM12-14 # 标准输出文件
#SBATCH -e logs/%j.treeM12-14 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree12.phy > tree12Model.outfile
ProteinModelSelection.pl tree13.phy > tree13Model.outfile
ProteinModelSelection.pl tree14.phy > tree14Model.outfile
echo end on $(date) # 结束时间
```

```tree_ModelSel_15-17.job
#!/bin/bash
#SBATCH -J treeM15-17 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeM15-17 # 标准输出文件
#SBATCH -e logs/%j.treeM15-17 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree15.phy > tree15Model.outfile
ProteinModelSelection.pl tree16.phy > tree16Model.outfile
ProteinModelSelection.pl tree17.phy > tree17Model.outfile
echo end on $(date) # 结束时间
```

```tree_ModelSel_18-20.job
#!/bin/bash
#SBATCH -J treeM18-20 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeM18-20 # 标准输出文件
#SBATCH -e logs/%j.treeM18-20 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
ProteinModelSelection.pl tree18.phy > tree18Model.outfile
ProteinModelSelection.pl tree19.phy > tree19Model.outfile
ProteinModelSelection.pl tree20.phy > tree20Model.outfile
echo end on $(date) # 结束时间
```


```bash
# (base) shichuang001@login01 一 4月 06 2026 14:02:36 ~/jlkang/Nocturnal_fish/trees
sbatch tree_ModelSel_3-5.job
# Submitted batch job 113821
sbatch tree_ModelSel_6-8.job
# Submitted batch job 113822
sbatch tree_ModelSel_9-11.job
# Submitted batch job 113823
sbatch tree_ModelSel_12-14.job
# Submitted batch job 113824
sbatch tree_ModelSel_15-17.job
# Submitted batch job 113825
sbatch tree_ModelSel_18-20.job
# Submitted batch job 113826
```

```tree_Cons_3-5.job
#!/bin/bash
#SBATCH -J treeC3-5 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeC3-5 # 标准输出文件
#SBATCH -e logs/%j.treeC3-5 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree3.phy -o Zebrafish -n tree3 -T 192 # 执行命令
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree4.phy -o Zebrafish -n tree4 -T 192
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree5.phy -o Zebrafish -n tree5 -T 192
echo end on $(date) # 结束时间
```

```tree_Cons_6-8.job
#!/bin/bash
#SBATCH -J treeC6-8 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeC6-8 # 标准输出文件
#SBATCH -e logs/%j.treeC6-8 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree6.phy -o Zebrafish -n tree6 -T 192 # 执行命令
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree7.phy -o Zebrafish -n tree7 -T 192
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree8.phy -o Zebrafish -n tree8 -T 192
echo end on $(date) # 结束时间
```

```tree_Cons_9-11.job
#!/bin/bash
#SBATCH -J treeC9-11 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeC9-11 # 标准输出文件
#SBATCH -e logs/%j.treeC9-11 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree9.phy -o Zebrafish -n tree9 -T 192 # 执行命令
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree10.phy -o Zebrafish -n tree10 -T 192
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree11.phy -o Zebrafish -n tree11 -T 192
echo end on $(date) # 结束时间
```

```tree_Cons_12-14.job
#!/bin/bash
#SBATCH -J treeC12-14 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeC12-14 # 标准输出文件
#SBATCH -e logs/%j.treeC12-14 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree12.phy -o Zebrafish -n tree12 -T 192 # 执行命令
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree13.phy -o Zebrafish -n tree13 -T 192
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree14.phy -o Zebrafish -n tree14 -T 192
echo end on $(date) # 结束时间
```

```tree_Cons_15-17.job
#!/bin/bash
#SBATCH -J treeC15-17 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeC15-17 # 标准输出文件
#SBATCH -e logs/%j.treeC15-17 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree15.phy -o Zebrafish -n tree15 -T 192 # 执行命令
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree16.phy -o Zebrafish -n tree16 -T 192
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree17.phy -o Zebrafish -n tree17 -T 192
echo end on $(date) # 结束时间
```

```tree_Cons_18-20.job
#!/bin/bash
#SBATCH -J treeC18-20 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.treeC18-20 # 标准输出文件
#SBATCH -e logs/%j.treeC18-20 # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree18.phy -o Zebrafish -n tree18 -T 192 # 执行命令
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree19.phy -o Zebrafish -n tree19 -T 192
raxmlHPC-SSE3 -f a -m PROTGAMMAJTTF -p 12345 -x 12345 -# 100 -s tree20.phy -o Zebrafish -n tree20 -T 192
echo end on $(date) # 结束时间
```

```bash
# (base) shichuang001@login01 一 4月 06 2026 19:54:11 ~/jlkang/Nocturnal_fish/trees
rm *phy_EVAL *phy_EVAL.out ST_tree*.phy_out RAxML_parsimonyTree.ST_tree*.phy
# (base) shichuang001@login01 一 4月 06 2026 19:59:11 ~/jlkang/Nocturnal_fish/trees
sbatch tree_Cons_3-5.job
# Submitted batch job 113838
sbatch tree_Cons_6-8.job
# Submitted batch job 113839
sbatch tree_Cons_9-11.job
# Submitted batch job 113840
sbatch tree_Cons_12-14.job
# Submitted batch job 113841
sbatch tree_Cons_15-17.job
# Submitted batch job 113842
sbatch tree_Cons_18-20.job
# Submitted batch job 113843

# prepare the input for paml
# h2076@h2076 Fri Apr 10 2026 16:37:24 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
mkdir paml_new
perl temp8.pl
```

```temp8.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

my (%hash1, %hash2);
my $name="Orth_";
my $info="orthologous_list_rep_2.txt";
open INFO, $info or die "can not open $info\n";
while (<INFO>) {
    chomp;
    my @a=split;
    next if /^Orth/;
    for (my $i = 1; $i < @a; $i++) {
        $hash1{$a[0]}.=$a[$i]."\t" unless $a[$i] eq "None";
    }
    $hash1{$a[0]}=~s/\s+$//;
}


my $tree="paml_rep_2_genes_treeInfo.txt";
open TREE, $tree or die "can not open $tree\n";
while (<TREE>) {
    chomp;
    my @a=split;
    $hash2{$a[-1]}++;
    my $file="paml_new/$name".$a[-1].".txt";
    open FILE, ">>$file" or die "can not create $file\n";
    if ($hash2{$a[-1]} == 1) {
        $a[1]=~s/\&\&/\t/g;
        print FILE "orth\t$a[1]\n";
        print FILE "$a[0]\t$hash1{$a[0]}\n";
    } else {
        print FILE "$a[0]\t$hash1{$a[0]}\n";
    }
}
```

```bash
# Create correlation file
# h2076@h2076 Fri Apr 10 2026 22:27:05 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new
perl temp1.pl
perl temp2.pl
perl temp3.pl
perl temp4.pl > final_paml_orth.txt

# install hyphy paml genewise
conda install hyphy -c bioconda
# activate hyphy
conda activate base

# h2076@h2076 Mon Apr 13 2026 22:07:35 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new/trees
perl temp5.pl

# h2076@h2076 Fri Apr 10 2026 22:27:05 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new
# hyphy detect the positive selected genes in nocturnal fish species
# Sublime edit the tree
# Replace: \)\:0\.\d+ => \)
# Replace: \:0\.\d+ => {Foreground}
nohup perl run_hyphy.pl final_paml_orth.txt Nocturnal > run_hyphy.process 2>&1 &
# h2076@h2076 Mon Apr 13 2026 23:45:16 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new/Hyphy
perl Extract_hyphy_Pvalue.pl
# no positive selection sites detected
# run Hyphy to detect the relax selection
nohup perl run_HyphyRelax.pl final_paml_orth.txt Nocturnal > run_HyphyRelax.process 2>&1 &
# h2076@h2076 Mon Apr 13 2026 23:46:22 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new
perl extract_relax.pl > Hyphy_relax.txt
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

my @txts=<*.txt>;
foreach my $txt (@txts) {
    my $tree;
    ($tree)=$txt=~/Orth_(.*?)\.txt/;
    print "$txt\t$tree\n";
    my $cor="correlation_".$tree.".txt";
    open COR, ">$cor" or die "can not create $cor\n";
    open TXT, $txt or die "can not open $txt\n";
    while (<TXT>) {
        chomp;
        my @a=split;
        if (/^orth/) {
            for (my $i = 1; $i < @a; $i++) {
                my $nuc=$a[$i]."_nuc.fasta";
                my $pep=$a[$i]."_pep.fasta";
                print COR "$a[$i]\t$pep\t$nuc\n";
            }
        }
    }
}
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

my @txts=<Orth_*.txt>;
foreach my $txt (@txts) {
    my ($tree, $dir);
    ($tree)=$txt=~/Orth_(.*?)\.txt/;
    ($dir)=$txt=~/(.*?)\.txt/;
    my $cor="correlation_".$tree.".txt";
    system("mkdir $dir");
    system("mv $txt $cor $dir/");
    system("cp prepare_input_paml.pl $dir/");
}
```

```temp3.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

for (my $i = 1; $i <= 20; $i++) {
    my $dir="Orth_tree"."$i";
    my $list=$dir.".txt";
    my $cor="correlation_tree".$i.".txt";
    chdir($dir);
    system("perl prepare_input_paml.pl --input $list --seq_dir ../../paml_input --cor_list $cor --output .");
    chdir("../");
}
```

```temp4.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

my %hash;
my $ano="../paml_rep_2_genes.txt";
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
    chomp;
    my @a=split /\t/;
    $hash{$a[0]}=$a[2]."\t".$a[3];
}

my @txts=<Orth_tree*/final_orth_input_paml.txt>;
foreach my $txt (@txts) {
    my $tree;
    ($tree)=$txt=~/(Orth_tree.*)\/final_orth_input_paml\.txt/;
    #print "$tree\n";
    my $info;
    open TXT, $txt or die "can not open $txt\n";
    while (<TXT>) {
        chomp;
        $info.=$_."\t".$hash{$_}."\t";
    }
    if ($info) {
        $info=~s/\s+$//;
        print "$tree\t$info\n";
    }
}
```

```run_hyphy.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# hyphy /lustre1/g/sbs_schunter/Kang/Ldim_revision/BUSTED-MH.bf --alignment final_alignment.fa --tree /lustre1/g/sbs_schunter/Kang/Ldim_revision/spe_hyphy.tre --branches Foreground
my $list=$ARGV[0]; # The list
my $spe =$ARGV[1]; # The species for positive selection analysis
my $outd="Hyphy";
unless (-d $outd) {
    mkdir $outd;
}

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
    chomp;
    my @a=split;
    my $orth=$a[0]."_".$a[1];
    my $alig="$a[0]/$a[1]/final_alignment.fa"; # the alignment
    my $outl="$outd/$orth"."_".$spe.".txt";
    my $tree="trees/$a[0]".".hyphy";
    my $cmd ="hyphy busted --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --branches Foreground > $outl";
#   print "$cmd\n";
    push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(30);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```run_HyphyRelax.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# hyphy /lustre1/g/sbs_schunter/Kang/Ldim_revision/BUSTED-MH.bf --alignment final_alignment.fa --tree /lustre1/g/sbs_schunter/Kang/Ldim_revision/spe_hyphy.tre --branches Foreground
my $list=$ARGV[0]; # The list
my $spe =$ARGV[1]; # The species for positive selection analysis
my $outd="Hyphy_relax";
unless (-d $outd) {
    mkdir $outd;
}

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
    chomp;
    my @a=split;
    my $orth=$a[0]."_".$a[1];
    my $alig="$a[0]/$a[1]/final_alignment.fa"; # the alignment
    my $outl="$outd/$orth"."_".$spe."_relax.txt";
    my $tree="trees/$a[0]".".hyphy";
    my $cmd ="hyphy relax --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --test Foreground > $outl";
#   print "$cmd\n";
    push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(30);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```extract_relax.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Array::Utils qw(:all);
use Parallel::ForkManager;

# final_alignment.fa.RELAX.json
my @relas=<Orth_tree*/OG*/final_alignment.fa.RELAX.json>;
foreach my $csu (@relas) {
     my $orth;
    ($orth)=$csu=~/(.*?)\/final_alignment\.fa\.RELAX\.json/;
    #print "$orth\n";
    my ($LRT, $p, $K);
    open CSU, $csu or die "can not open $csu\n";
    while (<CSU>) {
        chomp;
        if (/\"LRT\"\:(.*)\,/) {
            $LRT=$1;
        } elsif (/\"p-value\"\:(.*)\,/) {
                $p=$1;
        } elsif (/\"relaxation\s+or\s+intensification\s+parameter\"\:(.*)/) {
                $K=$1;
        } else {
                next;
        }
    }
    if ($LRT && $p && $K) {
        print "$orth\t$K\t$LRT\t$p\n";
    }
}
```

```extract_relax.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Array::Utils qw(:all);
use Parallel::ForkManager;

my %hash;
my $ortho="final_paml_orth.txt";
open ORTHO, $ortho or die "can not open $ortho\n";
while (<ORTHO>) {
    chomp;
    my @a=split /\t/;
    $hash{$a[1]}=$a[2]."\t".$a[3];
}

# final_alignment.fa.RELAX.json
my @relas=<Orth_tree*/OG*/final_alignment.fa.RELAX.json>;
foreach my $csu (@relas) {
    my $orth;
    my @a=split /\//, $csu;
    my $ano=$hash{$a[1]};

    ($orth)=$csu=~/(.*?)\/final_alignment\.fa\.RELAX\.json/;
    #print "$orth\n";
    my ($LRT, $p, $K);
    open CSU, $csu or die "can not open $csu\n";
    while (<CSU>) {
        chomp;
        if (/\"LRT\"\:(.*)\,/) {
            $LRT=$1;
        } elsif (/\"p-value\"\:(.*)\,/) {
                $p=$1;
        } elsif (/\"relaxation\s+or\s+intensification\s+parameter\"\:(.*)/) {
                $K=$1;
        } else {
                next;
        }
    }
    if ($LRT && $p && $K) {
        print "$orth\t$ano\t$K\t$LRT\t$p\n";
    }
}
```

## Detect the convergent sites
```bash
# 1. translate the DNA sequences to protein sequences
# h2076@h2076 Wed Apr 15 2026 07:00:42 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new
perl temp6.pl
# 2. Detect the convergent sites
perl Detect_Nons_all.pl > convergent_sites.txt

# Detect the specific genes among the orthologous genes
# h2076@h2076 Wed Apr 15 2026 07:56:52 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
grep 'CRY' all_swissprot_diamond_ano_final.txt|perl -alne '@a=split /\t/;my $nm;($nm)=$a[1]=~/sp\|.*\|(.*)\_.*/;$a[-1]=~s/\s+\[.*\]//;print "$nm\t$a[-1]"'|sort -u
```

```temp6.pl
#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_paml_orth.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
    chomp;
    my @a=split;
    my $orthdir=$a[0]."/".$a[1];
    chdir "$orthdir";
    my $cmd1="translateDna.pl -i final_alignment.fa > final_alignment_pep.fa";
    system($cmd1);
    chdir "/home/h2076/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_new/";
}
```

```Detect_Nons.pl
#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

# Detect_Nons.pl
my %seq; my $name;
my $orth=$ARGV[0];
my $fas="$orth/final_alignment_pep.fa";
open FAS, $fas or die "can not open $fas\n";
while (<FAS>) {
    chomp;
    if (/^>/) {
        s/\>//;
        $name=$_;
    } else {
        $seq{$name}.=$_;
    }
}

my %cleaner=(
    'Acura'=> 1,
    'Apoly'=> 1,
    'Daru'=> 1,
    'Pmol'=> 1,
    'Padel'=> 1,
    'Platyfish'=> 1,
    'Fugu'=> 1,
    'Medaka'=> 1,
    'Stickleback'=> 1,
    'Zebrafish'=> 1,
    );

# compare the nonsynonymous position pep sequences one by one
my %hash1;
my @poss;
my @nocls=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer);
my @cleas=qw(Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish);
my @aspes=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish);

&Build_pos_hash(\@aspes);

sub Build_pos_hash {
    my ($grp)=@_;
    my @grp=@{$grp};
    my $len;
    foreach my $spe (@grp) {
        if ($seq{$spe}) {
            my $seq=$seq{$spe};
            $len=length($seq);
            for (my $i = 0; $i < $len; $i++) {
                my $spepos=substr($seq,$i,1);
                $hash1{$spe}->{$i}=$spepos;
            }
        }
    }
    for (my $i = 0; $i < $len; $i++) {
        my (%hash2,%hash3);
        my $pos=$i;
        my $newp=$pos+1;
        my $info=$newp.":";
        foreach my $spe (@aspes) {
            if ($hash1{$spe}->{$pos}) {
                my $spepos=$hash1{$spe}->{$pos};
                $info.=$spe."($spepos);";
            }
        }

        my (@cleas_pos, @nocls_pos);
        foreach my $spe (@cleas) {
            if ($hash1{$spe}->{$pos}) {
                my $spepos=$hash1{$spe}->{$pos};
                $hash2{$spepos}++;
                push @cleas_pos, $spepos;
            }
        }
        foreach my $spe (@nocls) {
            if ($hash1{$spe}->{$pos}) {
                my $spepos=$hash1{$spe}->{$pos};
                $hash3{$spepos}++;
                push @nocls_pos, $spepos;
            }
        }
        my @isect = intersect(@cleas_pos, @nocls_pos);
        my $numb2=keys %hash2;
        my $numb3=keys %hash3;
        unless (@isect) {
            print "$orth\t$numb2\t$numb3\t$info\n" if $numb2==1 || $numb3==1;
        }
    }
}
```

```Detect_Nons_all.pl
#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_paml_orth.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
        chomp;
        my @a=split;
        my $name=$a[0]."/".$a[1];
        system("perl Detect_Nons.pl $name");
}
```








