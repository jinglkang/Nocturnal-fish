# Evolution of core clock genes (CCGs) and opsins in nocturnal fish
## Extract the orthogroups if any of zebrafish genes was annotated as target genes
### Target genes: opsins, core clock genes
```bash
# h2076@h2076 Fri Sep 18 2026 22:58:26 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
# 要求orthogroup里面必须有斑马鱼的序列，且该序列被注释为目标基因的名字，同时夜行鱼数量大于6，日行鱼数量大于10
perl Extract_target_orthogroups.pl > Target_Orthogroups.txt
# 包含所有38个物种序列的orthogroups已经被检查过了，所以此次关注小于38个物种的orthogroups
less Target_Orthogroups.txt|perl -alne 'next if /Orthogroup/i;my $nb=$F[3];print $F[4] if $nb<38' > Target_orthologous_list.txt
perl create_orth_rep.pl > Target_orthologous_list_rep.txt
mkdir Target_genes; cd Target_genes
# h2076@h2076 Fri Sep 18 2026 23:33:09 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/Target_genes
cp ../paml_input/*.fasta ./; cp ../paml_input/correlation.txt ./; cp ../Target_orthologous_list_rep.txt ./
cp ../paml_input/prepare_input_paml.pl ./
# perl prepare_input_paml.pl --input Target_orthologous_list_rep.txt --seq_dir . --cor_list correlation.txt --output .
nohup perl prepare_input_paml_parallel.pl Target_orthologous_list_rep.txt >prepare_input_paml.process 2>&1 & # 15 genes
wc -l final_orth_input_paml.txt 
# 剩下6个基因：OG0005031 (OPSD), OG0002759 (BHE40) , OG0005661 (NFIL3), OG0003231 (HLF), OG0014738 (OPSR2), OG0012368 (NFIL3)


# 使用Hyphy进行正选择和Relax选择压力分析，CASStools进行位点趋同分析
# 使用R只保留系统发育树的目标物种
# less final_alignment.fa|grep '>'|perl -alne 's/>//;print' > keep_species.txt
perl temp1.pl
# 添加前景枝标记（所有的夜行鱼类均作为前景枝），用于Hyphy分析
# h2076@h2076 Sat Sep 19 2026 12:16:41 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/Target_genes/OG0005031
vi hyphy.tre
# (Zebrafish,((Fugu,((Medaka,Platyfish),(Daru,((Acura,Apoly),(Pmol,Padel))))),(Tzosterophora{Foreground},(Tfucata{Foreground},(Rgracilis{Foreground},(((Snematoptera{Foreground},(Pfraenatus{Foreground},Pexostigma{Foreground})),((Fvariegata{Foreground},Acrassiceps{Foreground}),((Amelas{Foreground},Abrevicaudatus{Foreground}),((Nsavayensis{Foreground},Nviria{Foreground}),(Pmirifica{Foreground},Nfusca{Foreground}))))),((Fthermalis{Foreground},(Zleptacanthus{Foreground},Zviridiventer{Foreground})),((Cartus{Foreground},Cmacrodon{Foreground}),(((Onigrofasciatus{Foreground},Onovemfasciatus{Foreground}),(Ocookii{Foreground},Odoederleini{Foreground})),(Onotatus{Foreground},(Cquinquelineatus{Foreground},(Ocompressus{Foreground},(Ocyanosoma{Foreground},Oangustatus{Foreground})))))))))))));

nohup hyphy busted --alignment final_alignment.fa --tree hyphy.tre --multiple-hits Double+Triple --starting-points 5 --branches Foreground > hyphy_busted_results.txt 2>&1 &
nohup hyphy relax --alignment final_alignment.fa --tree hyphy.tre --multiple-hits Double+Triple --starting-points 5 --test Foreground > hyphy_relax_results.txt 2>&1 &

# Relax: OG0005661 (NFIL3), OG0012368 (NFIL3)

# 趋同位点
# 先翻译成蛋白序列
# h2076@h2076 Sat Sep 19 2026 15:11:21 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/Target_genes/OG0005031
translateDna.pl -i final_alignment.fa > final_alignment_pep.fa
less final_alignment_pep.fa|grep '>'|perl -alne 's/\>//;print'
vi config.tab
ct discovery -a final_alignment_pep.fa -t config.tab -o discovery.output --fmt fasta

# 整合成一个脚本
perl temp2.pl # 未发现这些基因存在patter1的convergence
```

```Create_newtree.R
library(ape)
tree <- read.tree("../spe.tre")
keep_species <- readLines("keep_species.txt")
subtree <- keep.tip(tree, keep_species)
write.tree(subtree, "subset_tree.nwk")
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(getcwd);

my $pwd=getcwd();
#print "$pwd\n";
my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
    chomp;
    my $orthid=$_;
    chdir($orthid);
    my $pwd1=getcwd();
    #print "$pwd1\n";
    my $cmd1="cp ../Create_newtree.R ./";
    system "$cmd1";
    my $fast="final_alignment.fa";
    open FAST, $fast or die "can not open $fast\n";
    my $keep="keep_species.txt";
    open KEEP, ">$keep" or die "can not create $keep\n";
    while (<FAST>) {
        chomp;
        if (/\>/) {
            s/\>//;
            print KEEP "$_\n";
        }
    }
    my $cmd2="Rscript Create_newtree.R";
    system($cmd2);
    chdir "$pwd";
    my $pwd2=getcwd();
    #print "$pwd2\n";
}
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(getcwd);

my %dius=&build_dius();
my $pwd=getcwd();
#print "$pwd\n";
my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
    chomp;
    my $orthid=$_;
    chdir($orthid);
    my $pwd1=getcwd();
    #print "$pwd1\n";
    my $cmd1="translateDna.pl -i final_alignment.fa > final_alignment_pep.fa";
    system "$cmd1";
    my $fast="final_alignment.fa";
    open FAST, $fast or die "can not open $fast\n";
    my $keep="config.tab";
    open KEEP, ">$keep" or die "can not create $keep\n";
    while (<FAST>) {
        chomp;
        if (/\>/) {
            s/\>//;
            my $tag;
            ($dius{$_})?($tag=0):($tag=1);
            print KEEP "$_\t$tag\n";
        }
    }
    my $cmd2="ct discovery -a final_alignment_pep.fa -t config.tab -o discovery.output --fmt fasta";
    system($cmd2);
    chdir "$pwd";
    my $pwd2=getcwd();
    #print "$pwd2\n";
}

sub build_dius {
    my %hash=(
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
    return(%hash);
}
```
