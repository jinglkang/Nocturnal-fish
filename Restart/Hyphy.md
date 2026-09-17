## estimate dN/dS for all branch based on the sequences by concatenated all single copy orthologous genes
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%hash1, %hash2); my @spes;
my $orthpep="../final_orth_input_paml.txt";
open ORTHPEP, $orthpep or die "can not open $orthpep\n";
while (<ORTHPEP>) {
	chomp;
	my $orthdir=$_;
	my $pep="../$orthdir/final_alignment.fa";
	my $spe;
	open PEP, $pep or die "can not open $pep\n";
	while (<PEP>) {
		chomp;
		if (/\>/) {
			s/\>//;
			$spe=$_;
			$hash1{$spe}++;
			push @spes, $spe if $hash1{$spe}==1;
		} else {
			$hash2{$spe}.=$_;
		}
	}
}

foreach my $spe (@spes) {
	print ">$spe\n$hash2{$spe}\n";
}
```

```bash
# perl codeml.pl --input temp/$temp --model free-ratio --dir . --tree spe.tre --icode 0 --omega 1.2
# jlkang@hnu2024 Thu Dec 04 2025 20:18:33 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
mkdir free_ratio_conca
cp OG0013142/free-ratio.ctr free_ratio_conca/
cp OG0013142/spe.tre free_ratio_conca/
# jlkang@hnu2024 Thu Dec 04 2025 20:41:07 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/free_ratio_conca
perl temp1.pl > final_alignment_conca.fa
fasta2phy.pl final_alignment_conca.fa; mv final_alignment_conca.fa.phy final_alignment_conca.phy
vi free-ratio.ctr
nohup codeml free-ratio.ctr > free_ratio.process 2>&1 &
```

## Detect the positive selection signals in nocturnal fishes
```run_hyphy.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# hyphy /lustre1/g/sbs_schunter/Kang/Ldim_revision/BUSTED-MH.bf --alignment final_alignment.fa --tree /lustre1/g/sbs_schunter/Kang/Ldim_revision/spe_hyphy.tre --branches Foreground
my $list=$ARGV[0]; # The list
my $tree=$ARGV[1]; # The tree
my $spe =$ARGV[2]; # The species for positive selection analysis
my $outd="Hyphy";
unless (-d $outd) {
        mkdir $outd;
}

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        my @a=split;
        my $orth=$a[0];
        my $alig="$orth/final_alignment.fa"; # the alignment
        my $outl="$outd/$orth"."_".$spe.".txt";
        my $cmd ="hyphy busted --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --branches Foreground > $outl";
#       print "$cmd\n";
        push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(110);
foreach my $cmd (@cmds) {
        $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```bash
# jlkang@hnu2024 Fri Nov 14 2025 11:10:08 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
vi spe_hyphy.tre
# (((((Snematoptera{Foreground},(Pexostigma{Foreground},Pfraenatus{Foreground})),(((Tfucata{Foreground},(Tzosterophora{Foreground},((Zebrafish,Stickleback),(Fugu,((Platyfish,Medaka),(((Padel,Pmol),(Apoly,Acura)),Daru)))))),Rgracilis{Foreground}),(((Zviridiventer{Foreground},Zleptacanthus{Foreground}),Fthermalis{Foreground}),((((Odoederleini{Foreground},Ocookii{Foreground}),(Onovemfasciatus{Foreground},Onigrofasciatus{Foreground})),(Onotatus{Foreground},((Ocompressus{Foreground},(Oangustatus{Foreground},Ocyanosoma{Foreground})),Cquinquelineatus{Foreground}))),(Cmacrodon{Foreground},Cartus{Foreground}))))),(Acrassiceps{Foreground},Fvariegata{Foreground})),((Nfusca{Foreground},Pmirifica{Foreground}),(Nviria{Foreground},Nsavayensis{Foreground}))),Amelas{Foreground},Abrevicaudatus{Foreground});
# run hyphy with BUSTED-MH method
# final_orth_input_paml.txt
conda activate hyphy_env
# hyphy busted --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --branches Foreground
nohup perl run_hyphy.pl final_orth_input_paml.txt spe_hyphy.tre Nocturnal > run_hyphy.process 2>&1 &
# [1] 90730

# (base) jlkang@hnu2024 Thu Nov 20 2025 13:34:57 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/Hyphy
perl Extract_hyphy_Pvalue.pl > Hyphy_pValue.txt
R
# p_apoly<-read.table(file="Hyphy_pValue.txt")
# p_apoly$fdr<- p.adjust(p_apoly$V2,method="fdr",length(p_apoly$V2))
# write.table(p_apoly, file="Hyphy_pValue_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
less Hyphy_pValue_fdr.txt | perl -alne 'print if $F[2]<=0.05' > Hyphy_PSGs.txt
# jlkang@hnu2024 Fri Dec 05 2025 10:12:01 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl anno_orth.pl > Hyphy_PSGs_anno.txt
```

## Detect the strength of selection
```run_HyphyRelax.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# hyphy /lustre1/g/sbs_schunter/Kang/Ldim_revision/BUSTED-MH.bf --alignment final_alignment.fa --tree /lustre1/g/sbs_schunter/Kang/Ldim_revision/spe_hyphy.tre --test Foreground
my $list=$ARGV[0]; # The list
my $tree=$ARGV[1]; # The tree
my $outd="Hyphy_relax";
unless (-d $outd) {
        mkdir $outd;
}

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        my @a=split;
        my $orth=$a[0];
        my $alig="$orth/final_alignment.fa"; # the alignment
        my $outl="$outd/$orth"."_relax.txt";
        my $cmd ="hyphy relax --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --test Foreground > $outl";
#       print "$cmd\n";
        push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(110);
foreach my $cmd (@cmds) {
        $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```bash
# Hyphy relax
conda activate hyphy_env
# hyphy relax --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --branches Foreground
nohup perl run_HyphyRelax.pl final_orth_input_paml.txt spe_hyphy.tre > run_HyphyRelax.process 2>&1 &
# [1] 502469
# jlkang@hnu2024 Thu Dec 04 2025 17:25:04 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0013632
less final_alignment.fa.RELAX.json|grep -i 'relaxation or intensification parameter'
#   "relaxation or intensification parameter":1.658051047455258
# "test results":{
#   "LRT":12.05207283040363,
#   "p-value":0.0005173482067816204,
#   "relaxation or intensification parameter":1.658051047455258
#  },
# jlkang@hnu2024 Mon Dec 08 2025 11:58:44 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl extract_relax.pl > Hyphy_relax_genes.txt
# FDR
R
# p_apoly<-read.table(file="Hyphy_relax_genes.txt")
# p_apoly$fdr<- p.adjust(p_apoly$V4,method="fdr",length(p_apoly$V4))
# write.table(p_apoly, file="Hyphy_relax_genes_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
less Hyphy_relax_genes_fdr.txt|perl -alne 'print if $F[-1]<=0.05' > Hyphy_relax_genes_fdr_sig.txt
perl anno_orth.pl > Hyphy_relax_genes_fdr_sig_ano.txt
less Hyphy_relax_genes_fdr_sig.txt|perl -alne 'print "$F[0]" if $F[1]>1' > Intensified_genes.txt
less Hyphy_relax_genes_fdr_sig.txt|perl -alne 'print "$F[0]" if $F[1]<1' > Relaxed_genes.txt
```

```extract_relax.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Array::Utils qw(:all);
use Parallel::ForkManager;

# final_alignment.fa.RELAX.json
my @relas=<paml_input/OG*/final_alignment.fa.RELAX.json>;
foreach my $csu (@relas) {
        my $orth;
    ($orth)=$csu=~/paml_input\/(OG.*)\/final_alignment\.fa\.RELAX\.json/;
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

### Install csubst
```foreground.txt
1	Tfucata
1	Tzosterophora

2	Rgracilis

3	Fvariegata

4	Acrassiceps

5	Abrevicaudatus
5	Amelas
5	Nsavayensis
5	Nviria
5	Pmirifica
5	Nfusca
5	Snematoptera

6	Pfraenatus
6	Pexostigma

7	Fthermalis
7	Zviridiventer
7	Zleptacanthus

8	Cquinquelineatus
8	Cmacrodon
8	Cartus

9	Odoederleini
9	Ocookii
9	Onovemfasciatus
9	Onigrofasciatus
9	Onotatus
9	Ocompressus
9	Oangustatus
9	Ocyanosoma
```

```bash
# 源代码安装才能成功
# 需要安装iqtree2（但是不能确定是不是这个原因）
# 安装csubst分析convergence
下载安装包（https://github.com/kfuku52/csubst/releases/tag/v1.4.20）
解压
pip install csubst-1.4.20/
# csubst analyze --alignment_file alignment.fa --rooted_tree_file tree.nwk --foreground foreground.txt
# jlkang@hnu2024 Mon Dec 01 2025 22:41:11 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0000385
csubst analyze --alignment_file final_alignment_pep.fa --rooted_tree_file spe.tre --foreground foreground.txt
cp ../RAxML_bestTree.conca_sigpep_zebrafish_root ./
pip install git+https://github.com/kfuku52/nwkit
nwkit drop --infile RAxML_bestTree.conca_sigpep_zebrafish_root --target intnode --name yes --outfile tree.nwk
# use the nucleotide sequence fasta file
csubst analyze --alignment_file final_alignment.fa --rooted_tree_file RAxML_bestTree.conca_sigpep_zebrafish_root --foreground foreground.txt
# jlkang@hnu2024 Wed Dec 03 2025 12:03:56 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup perl run_csubst.pl final_orth_input_paml.txt > run_csubst.process 2>&1 &
# [1] 194770
less 1.txt |perl -alne 's/,//g;my @a=split;foreach my $i(@a){$j++;print "$j\t$i" if $i>=10}'
less final_alignment_pep.fa|perl -alne 'if (/>/){s/>//;$nm=$_;}else{my $id=substr($_,337,1);print "$nm\t$id"}'
for i in OG*/csubst_cb_2.tsv;do less ${i}|perl -alne 'if ($F[7] eq "Y" && $F[47]>=2 && $F[12]>=5){print;last}';done|wc -l
# 102 there are 120 genes with convergence in nocturnal fish species

# jlkang@hnu2024 Fri Dec 05 2025 11:16:00 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl extract_csubst.pl > csubst_genes.txt
perl anno_orth.pl > csubst_genes_anno.txt
# kangjingliang@KangdeMacBook-Pro-2 五 12 05 2025 11:32:08 ~/Documents/2025/Nocturnal_fish
less csubst_genes_anno.txt|perl -alne 'my @a=split /\t/;print "$a[0]\t$a[1]\t$a[2]"'|sort -u > csubst_genes_anno_uniq.txt

# make sure that the convergent site were not same with any amino acid at the same site within the diurnal fish species
# jlkang@hnu2024 Sun Dec 07 2025 17:22:25 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl csubst_convergent.pl > csubst_ccs_genes.txt
perl anno_orth.pl > csubst_ccs_genes_anno.txt
less csubst_ccs_genes_anno.txt|cut -f 1,2,3|sort -u > csubst_ccs_genes_anno_uniq.txt
less csubst_genes_anno.txt|perl -alne 'print $F[0]' > csubst_ccs_genes_id.txt
```

```run_csubst.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# csubst analyze --alignment_file final_alignment.fa --rooted_tree_file RAxML_bestTree.conca_sigpep_zebrafish_root --foreground foreground.txt
my $list=$ARGV[0]; # The list
my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        my @a=split;
        my $orth=$a[0];
        my $alig="final_alignment.fa";
        my $tree="RAxML_bestTree.conca_sigpep_zebrafish_root";
        my $fore="foreground.txt";
        system("cp $tree $orth/");
        system("cp $fore $orth/");
        my $cmd ="cd $orth; csubst analyze --alignment_file $alig --rooted_tree_file $tree --foreground $fore";
#       print "$cmd\n";
        push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(110);
foreach my $cmd (@cmds) {
        $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```extract_csubst.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# csubst_cb_2.tsv
my @csus=<paml_input/OG*/csubst_cb_2.tsv>;
foreach my $csu (@csus) {
        my $orth;
        ($orth)=$csu=~/paml_input\/(OG.*)\/csubst_cb_2\.tsv/;
        #print "$orth\n";
        open CSU, $csu or die "can not open $csu\n";
        while (<CSU>) {
                chomp;
                my @a=split /\t/;
                if ($a[7] eq "Y" && $a[47]>=2 && $a[12]>=5) {
                        print "$orth\t$a[0]\t$a[1]\t$a[12]\t$a[47]\n";
                }
        }
}
```

```csubst_convergent.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Array::Utils qw(:all);
use Parallel::ForkManager;

my @dius=qw(Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish);
# putative convergent branches
my $PuCo="csubst_genes.txt";
open PUCO, $PuCo or die "can not open $PuCo\n";
while (<PUCO>) {
        my @a=split;
        my ($orth, $bran1, $bran2)=($a[0], $a[1], $a[2]);
        my $subs="paml_input/$orth/csubst_b.tsv";
        my (@cov1, @cov2);
        my ($spe1, $spe2);
        open SUBS, $subs or die "can not open $subs\n";
        while (<SUBS>) {
                chomp;
                my @a=split;
                if ($a[1] eq $bran1) {
                        $spe1=$a[0];
                        @cov1=split /\,/, $a[4];
                } elsif ($a[1] eq $bran2) {
                        $spe2=$a[0];
                        @cov2=split /\,/, $a[4];
                } else {
                        next;
                }
        }

        my $Fcov;

        my @isecs=intersect(@cov1, @cov2);
        foreach my $isec (@isecs) {
                my ($posi, $subN)=$isec=~/\D(\d+)(\D)/;
                #print "$orth\t$spe1\t$spe2\t$posi\t$subN\n";
                my $fas="paml_input/$orth/final_alignment_pep.fa";
                my %hash=&build_hash($fas);
                my @amins;
                foreach my $diu (@dius) {
                        my $seq=$hash{$diu};
                        my $id =$posi-1;
                        my $amin=substr($seq,$id,1);
                        push @amins, $amin;
                }
                my $i;
                foreach my $amin (@amins) {
                        if ($amin eq $subN) {
                                $i++;
                                if ($i >= 1) {
                                        last;
                                }
                        }
                }
                unless ($i) {
                        $Fcov.=$isec.",";
                }
        }
        if ($Fcov) {
                $Fcov=~s/\,$//;
                print "$orth\t$spe1\t$spe2\t$Fcov\n";
        }
}

sub build_hash {
        my ($fas)=@_;
        my %seqs; my $spe;
        open FAS, $fas or die "can not open $fas\n";
        while (<FAS>) {
                chomp;
                if (/>/) {
                        s/>//;
                        $spe=$_;
                } else {
                        $seqs{$spe}.=$_;
                }
        }
        return %seqs;
}
```

### Install PCOC to detect the convergent evolution
```bash
docker pull carinerey/pcoc
# run in the working directory
alias CMD_PCOC_DOCKER="docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD -e CWD=$PWD carinerey/pcoc"
# jlkang@hnu2024 Wed Nov 26 2025 15:11:30 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0000385
CMD_PCOC_DOCKER pcoc_num_tree.py -t spe.tre -o num_tree.pdf
scenario="0,1,2,5,6,28,30,31,33,35,36,38,39,42,43,44,45,48,52,53,59,60,63,64,66,67,71,72"
CMD_PCOC_DOCKER pcoc_det.py -t tree.tre -aa ali.fa -o output_pcoc_det -m $scenario -f 0.8
CMD_PCOC_DOCKER pcoc_det.py -t tree.nw -aa ali.fa -o output_pcoc_det -m $scenario --plot_complete_ali --plot
```

## dN/dS comparison
```bash
# plot according the dN/dS values (from the great to small)
# kangjingliang@KangdeMacBook-Pro-2 三 12 10 2025 21:07:15 ~/Documents/2025/Nocturnal_fish/paml_FreeRatio
perl All_dNdS_order.pl > free_ratio_result_3.txt # plot
# dNdS_median.csv: plot the median value

# PSGs
# plot the the dN/dS values of PSGs according the dN/dS values (from the great to small)
# kangjingliang@KangdeMacBook-Pro-2 三 12 10 2025 21:10:23 ~/Documents/2025/Nocturnal_fish/paml_FreeRatio
perl PSGs_dNdS.pl > free_ratio_result_3_PSGs.txt
perl PSGs_dNdS_order.pl > free_ratio_result_3_PSGs_plot.txt # plot
# dNdS_median_PSGs.csv: plot the median value

# CCS genes
# kangjingliang@KangdeMacBook-Pro-2 三 12 10 2025 23:58:32 ~/Documents/2025/Nocturnal_fish/paml_FreeRatio
perl CCSgenes_dNdS.pl > free_ratio_result_3_CCSgenes.txt
perl CCSgenes_dNdS_order.pl > free_ratio_result_3_CCSgenes_plot.txt # plot
```

## Phototransduction genes
```bash
# Extract the seuqences of five diurnal ref fish species
# Zebrafish_gene_ensembl.txt; Stickleback_gene_ensembl.txt; Platyfish_gene_ensembl.txt; Medaka_gene_ensembl.txt; Fugu_gene_ensembl.txt
# kangjingliang@KangdeMBP-2 五 12 12 2025 19:27:55 ~/Documents/2025/Nocturnal_fish/Specific_genes
# Phototransduction_genes.txt: the name of target genes; 
# Phototransduction: the directory of output sequences files

# kangjingliang@KangdeMBP-2 五 12 12 2025 20:24:52 ~/Documents/2025/Nocturnal_fish/Specific_genes
perl extract_target_genes.pl Phototransduction_genes.txt Phototransduction

# download all pep sequences of phototransduction genes 
# kangjingliang@KangdeMBP-2 五 12 12 2025 20:26:20 ~/Documents/2025/Nocturnal_fish/Specific_genes/Phototransduction
cat *.fas > Phototransduction.fasta
# jlkang@hnu2024 Sat Dec 13 2025 18:19:49 ~/Nocturnal_fish/Orthologous/pep
mkdir Phototransduction
mv Zebrafish.fas Stickleback.fas Platyfish.fas Medaka.fas Fugu.fas Phototransduction/
cat *.fas > all_test.fas; mv Phototransduction/*.fas ./; mv all_test.fas Phototransduction/
# kangjingliang@KangdeMBP-2 六 12 13 2025 18:26:33 ~/Documents/2025/Nocturnal_fish/Specific_genes/Phototransduction
scp Phototransduction.fasta jlkang@10.33.247.14:~/Nocturnal_fish/Orthologous/pep/Phototransduction/

# blast all sequences of the test species to Phototransduction.fasta
# diamond
# jlkang@hnu2024 Sat Dec 13 2025 18:25:38 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
diamond makedb --in Phototransduction.fasta -d Phototransduction
diamond blastp -q all_test.fas -e 1e-5 --sensitive -k 1 -d Phototransduction --out Phototransduction_blastp.result

# prepare the fasta file per gene for the input of orthofinder
# jlkang@hnu2024 Sun Dec 14 2025 15:59:37 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl prepare_seq_orthofinder.pl Phototransduction_blastp.result Phototransduction.fasta
# run orthofinder for each gene
for i in *_gene;do cd ${i};orthofinder -f ./;cd ../;done
# select the orthogroups with 80% nocturnal (23) and diurnal (8) fish species transcripts
perl select_qualified_orth.pl > Qualified_orth.txt
less Qualified_orth.txt|cut -f 1|sort -u|wc -l # 30: no duplicated genes

# obtain the orthologous_list_rep.txt (only keep the representative transcript per species)

# grep 'Platyfish' grk1a_gene/OrthoFinder/Results_Dec14/Orthogroups/Orthogroups.tsv

# based on Qualified_orth.txt, create "orthologous_list_rep.txt" for each gene
# jlkang@hnu2024 Tue Dec 16 2025 14:51:46 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl create_orth_rep_list.pl Phototransduction_blastp.result Phototransduction.fasta

# prepare the input sequences for evolutionary analyse
# jlkang@hnu2024 Tue Dec 16 2025 20:27:13 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl create_orth_seq.pl Phototransduction_blastp.result # all sequences are saved in "sequences/"

# create the "correlation.txt" in "*_gene/OrthoFinder/Results_Dec14/Orthogroups"
# jlkang@hnu2024 Tue Dec 16 2025 21:37:56 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl create_correlation.pl

# jlkang@hnu2024 Tue Dec 16 2025 23:41:26 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl Recode.pl > Qualified_orth_recode.txt

# picked the required species for phylogenetic tree construction
# based on the previous orthologous genes protein alignment
# jlkang@hnu2024 Wed Dec 17 2025 16:37:49 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
cp ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/conca_sigpep.fa ./
fasta2phy.pl conca_sigpep.fa # obtain "conca_sigpep.fa.phy"
nohup perl Build_phylogeny.pl > tree.reports 2>&1 & # all the phylogeny will be constructed in "Phylogeny/"
# [1] 2916856
# run the phylogeny "Orth_36spe_10code.phy" which didn't run for the last time
# jlkang@hnu2024 Sat Dec 27 2025 13:28:22 ~/Nocturnal_fish/Orthologous/pep/Phototransduction/Phylogeny
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s Orth_36spe_10code.phy -n Orth_36spe_10code -T 190 > tree.reports 2>&1 &
# [1] 5477

# prepare the input for evolutionary analysis
# jlkang@hnu2024 Mon Dec 22 2025 00:51:39 ~/Nocturnal_fish/Orthologous/pep/Phototransduction
perl prepare_evo.pl

# create the corresponding species tree for each phototransduction gene
# jlkang@hnu2024 Mon Dec 22 2025 01:51:30 ~/Nocturnal_fish/Orthologous/pep/Phototransduction/Phylogeny
perl create_evo_tre.pl
# put the tree into the corresponding folder

```

```prepare_seq_orthofinder.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%hash, %hash1, %hash2, %hash3);
my (@spes, @genes);

my $blast=$ARGV[0];
open BLAST, $blast or die "can not open $blast\n";
while (<BLAST>) {
        chomp;
        my @a=split /\t/;
        my ($test_spe, $ref_spe, $gene);
        ($test_spe)=$a[0]=~/(.*?)\_.*/;
        ($ref_spe, $gene)=$a[1]=~/(.*)\_.*\_(.*)/;

        $hash1{$a[0]}++;
        push @{$hash{$gene}->{$test_spe}}, $a[0] if $hash1{$a[0]}==1;

        $hash1{$a[1]}++;
        push @{$hash{$gene}->{$ref_spe}}, $a[1] if $hash1{$a[1]}==1;

        $hash2{$test_spe}++;
        push @spes, $test_spe if $hash2{$test_spe}==1;

        $hash2{$ref_spe}++;
        push @spes, $ref_spe if $hash2{$ref_spe}==1;

        $hash2{$gene}++;
        push @genes, $gene if $hash2{$gene}==1;
}

my $test_id;
my $testfas="all_test.fas";
open TEST, $testfas or die "can not open $testfas\n";
while (<TEST>) {
        chomp;
        if (/\>/) {
                s/\>//;
                $test_id=$_;
        } else {
                $hash3{$test_id}.=$_;
        }
}

my $ref_id;
my $reffas=$ARGV[1];
open REF, $reffas or die "can not open $reffas\n";
while (<REF>) {
        chomp;
        if (/\>/) {
                s/\>//;
                $ref_id=$_;
        } else {
                $hash3{$ref_id}.=$_;
        }
}

foreach my $gene (@genes) {
        my $dir=$gene."_gene";
        system("mkdir $dir");
        foreach my $spe (@spes) {
                if ($hash{$gene}->{$spe}) {
                        my $testfa="$dir/$spe.fas";
                        open TESTFA, ">>$testfa" or die "can not open $testfa\n";
                        my @a=@{$hash{$gene}->{$spe}};
                        foreach my $a (@a) {
                                my $seq=$hash3{$a};
                                print TESTFA ">$a\n$seq\n";
                        }
                }
        }
}
```

```select_qualified_orth.pl
#!/usr/bin/perl
use strict;
use warnings;

my @nocs=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer);
my @dius=qw(Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish);

my %nocspe;
foreach my $noc (@nocs) {
        $nocspe{$noc}++;
}

my @counts=<*_gene/OrthoFinder/Results_Dec14/Orthogroups/Orthogroups.GeneCount.tsv>;
foreach my $count (@counts) {
        my ($gene)=$count=~/(.*)_gene\/OrthoFinder\/Results_Dec14\/Orthogroups\/Orthogroups\.GeneCount\.tsv/;
        &stab_nb($count, $gene);
}

sub stab_nb {
        my ($count, $gene)=@_;
        my @heads;
        open COUNT, $count or die "can not open $count\n";
        while (<COUNT>) {
                chomp;
                my @a=split;
                my ($nb_noc, $nb_diu, $info);
                if (/^Orthogroup/) {
                        @heads=@a;
                } else {
                        for (my $i = 1; $i < @a-1; $i++) {
                                if ($a[$i]>=1) {
                                        $info.=$heads[$i]."&&";
                                        $nocspe{$heads[$i]}?($nb_noc++):($nb_diu++);
                                        #print "$nb_noc\t$nb_diu\n";
                                }
                        }
                } #  && $nb_noc >= 19 && $nb_diu >= 5
                if ($nb_noc && $nb_diu && $nb_noc >= 23 && $nb_diu >= 8) {
                        $info=~s/&&$//;
                        print "$gene\t$a[0]\t$nb_noc\t$nb_diu\t$info\n";
                }
        }
}
```

```create_orth_rep_list.pl
#!/usr/bin/perl
use strict;
use warnings;

my $blast=$ARGV[0];
my %hash1;
open BLAST, $blast or die "can not open $blast\n";
while (<BLAST>) {
    chomp;
    my @a=split /\t/;
        $hash1{$a[0]}=$a[-2]; # %hash1: eValue    
}

my @tests=qw(Acura Apoly Daru Pmol Padel Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer);
my @refs=qw(Platyfish Fugu Medaka Stickleback Zebrafish);

my %testspe;
foreach my $test (@tests) {
        $testspe{$test}++;
}

my %refspe;
foreach my $ref (@refs) {
        $refspe{$ref}++;
}

my $whole; my %seqs;
my $fasta=$ARGV[1]; # the fasta file of target genes
open FASTA, $fasta or die "can not open $fasta\n";
while (<FASTA>) {
        chomp;
        if (/>/) {
                s/>//;
                $whole=$_;
        } else {
                $seqs{$whole}.=$_;
        }
}

my $targ="Qualified_orth.txt";
open TARG, $targ or die "can not open $targ\n";
while (<TARG>) {
        chomp;
        my %hash2;
        my @headers;
        my @a=split /\t/;
        $hash2{$a[1]}++; # %hash2 for orth id of this orthogroup
        my $gene_nm=$a[0];

        my $tsv="$a[0]"."_gene/OrthoFinder/Results_Dec14/Orthogroups/Orthogroups.tsv";
        #my $header="Gene_nm\tOrthogroup\t";
        my $header="Orthogroup\t";
        my ($info, $number);
        open TSV, $tsv or die "can not open $tsv\n";
        while (<TSV>) {
                chomp;
                s/\r//g;
                s/\,//g;
                my @b=split /\t/;
                if (/^Orthogroup/) {
                        @headers=@b;
                        $number=@headers;
                } else {
                        if ($hash2{$b[0]}) {
                                #$info=$gene_nm."\t".$b[0]."\t";
                                $info=$b[0]."\t";
                                for (my $i = 1; $i < $number; $i++) {
                                        print "$gene_nm\t$headers[$i]\t$b[$i]\n" if $i==$number-1;
                                        if ($b[$i]) {
                                                my @c=split /\s+/, $b[$i];
                                                my $nb=@c;
                                                if (@c==1) {
                                                        $header.=$headers[$i]."\t";
                                                        $info.=$b[$i]."\t";                                                     
                                                        #print "$gene_nm\t$b[$i]\t$nb\n";
                                                } elsif (@c > 1 && $testspe{$headers[$i]}) {
                                                        #&test_spe_element($i,\@c);
                                                        my %hash3;
                                                        foreach my $c (@c) {
                                                                #print "$gene_nm\t$c\t$nb\n";
                                                                if ($hash3{$headers[$i]}) {
                                                                        my $old_evalue=$hash3{$headers[$i]}->{'EVE'};
                                                                        my $new_evalue=$hash1{$c};
                                                                        unless ($new_evalue || $old_evalue) {
                                                                                die "$gene_nm\t$headers[$i]\t$tsv\tTESTFISH\n";
                                                                        }
                                                                        if ($new_evalue < $old_evalue) {
                                                                                my $evalue=$hash1{$c};
                                                                                $hash3{$headers[$i]}={
                                                                                        'EVE' => $evalue,
                                                                                        'TRA' => $c
                                                                                };
                                                                        }
                                                                } else {                                        
                                                                        my $evalue=$hash1{$c};
                                                                        $hash3{$headers[$i]}={
                                                                                'EVE' => $evalue,
                                                                                'TRA' => $c
                                                                        };
                                                                }
                                                        }
                                                        $header.=$headers[$i]."\t";
                                                        $info.=$hash3{$headers[$i]}->{'TRA'}."\t";
                                                } elsif (@c > 1 && $refspe{$headers[$i]}) {
                                                        my %hash4;
                                                        foreach my $c (@c) {
                                                                if ($hash4{$headers[$i]}) {
                                                                        my $oldLen=$hash4{$headers[$i]}->{'LEN'};
                                                                        #print "$gene_nm\t$oldLen\t$headers[$i]\n";
                                                                        my $newLen=length($seqs{$c});
                                                                        #print "$gene_nm\t$newLen\t$headers[$i]\n";
                                                                        unless ($oldLen || $newLen) {
                                                                                #print "$headers[$i]\n";
                                                                                #print "REFFISH\t$headers[$i]\n";
                                                                                die "$gene_nm\t$tsv\tREFFISH\t$headers[$i]\n";
                                                                        }
                                                                        if ($newLen > $oldLen) {
                                                                                $hash4{$headers[$i]}={
                                                                                        'LEN' => $newLen,
                                                                                        'TRA' => $c
                                                                                };
                                                                        }
                                                                } else {
                                                                        my $len=length($seqs{$c});
                                                                        $hash4{$headers[$i]}={
                                                                                'LEN' => $len,
                                                                                'TRA' => $c
                                                                        };
                                                                }
                                                        }
                                                        $header.=$headers[$i]."\t";
                                                        $info.=$hash4{$headers[$i]}->{'TRA'}."\t";
                                                }       
                                        }
                                }
                                $header=~s/\s+$//;
                                $info=~s/\s+$//;
                                my $rep="$a[0]"."_gene/OrthoFinder/Results_Dec14/Orthogroups/orthologous_list_rep.txt";
                                open REP, ">>$rep" or die "can not create $rep\n";
                                print REP "$header\n$info\n";
                        }
                }
        }
}
```

```create_orth_seq.pl
#!/usr/bin/perl
use strict;
use warnings;

system("cat ~/Nocturnal_fish/Orthologous/pep/*.fas > all_pep.fas");
system("cat ~/Nocturnal_fish/Orthologous/nuc/*.fas > all_nuc.fas");

my (%pep, %nuc);
%pep=&build_hash("all_pep.fas");
%nuc=&build_hash("all_nuc.fas");

system("mkdir sequences");
my $blast=$ARGV[0];
my %hash1;
open BLAST, $blast or die "can not open $blast\n";
while (<BLAST>) {
    chomp;
    my @a=split /\t/;
        my ($testspe, $refspe);

        ($testspe)=$a[0]=~/(.*)\_.*/;
        my $peps="sequences/".$testspe."_pep.fasta";
        open PEPS, ">>$peps" or die "can not create $peps\n";
        print PEPS ">$a[0]\n$pep{$a[0]}\n";
        my $nucs="sequences/".$testspe."_nuc.fasta";
        open NUCS, ">>$nucs" or die "can not create $nucs\n";
        print NUCS ">$a[0]\n$nuc{$a[0]}\n";


        ($refspe)=$a[1]=~/(.*)\_.*\_.*/;
        (my $id) =$a[1]=~/(.*)\_.*/;
        print "$refspe\t$id\n";
        $hash1{$a[1]}++;
        if ($hash1{$a[1]}==1) {
                my $refpeps="sequences/".$refspe."_pep.fasta";
                open REFPEPS, ">>$refpeps" or die "can not create $refpeps\n";
                print REFPEPS ">$a[1]\n$pep{$id}\n";
                my $refnucs="sequences/".$refspe."_nuc.fasta";
                open REFNUCS, ">>$refnucs" or die "can not create $refnucs\n";
                print REFNUCS ">$a[1]\n$nuc{$id}\n";
        }
}

sub build_hash {
        my ($fas)=@_;
        my %hash;
        my $name;
        open FAS, $fas or die "can not open $fas\n";
        while (<FAS>) {
                chomp;
                if (/>/) {
                        s/>//;
                        $name=$_;
                } else {
                        $hash{$name}.=$_;
                }
        }
        return %hash;
}
```

```create_correlation.pl
#!/usr/bin/perl
use strict;
use warnings;

my $qua="Qualified_orth.txt";
open QUA, $qua or die "can not open $qua\n";
while (<QUA>) {
        chomp;
        my @a=split; my $j;
        my $orth=$a[0]."_gene/OrthoFinder/Results_Dec14/Orthogroups/orthologous_list_rep.txt";
        my $core=$a[0]."_gene/OrthoFinder/Results_Dec14/Orthogroups/correlation.txt";
        open CORE, ">$core" or die "can not create $core\n";
        open ORTH, $orth or die "can not open $orth\n";
        while (<ORTH>) {
                chomp;
                my $j++;
                my @b=split;
                for (my $i = 1; $i < @b; $i++) {
                        print CORE "$b[$i]\t$b[$i]_pep.fasta\t$b[$i]_nuc.fasta\n";
                }
                last if $j==1;
        }
}
```

```Recode.pl
#!/usr/bin/perl
use strict;
use warnings;

my $qua="Qualified_orth.txt";
my (%hash1, %hash2);
open QUA, $qua or die "can not open $qua\n";
while (<QUA>) {
        chomp;
        my @a=split;
        $hash1{$a[-1]}++;
}

my $i;
foreach my $key (sort keys %hash1) {
        $i++;
        $hash2{$key}=$i;
}

open QUA, $qua or die "can not open $qua\n";
while (<QUA>) {
        chomp;
        my @a=split;
        print "$_\t$hash2{$a[-1]}\n";
}
```

```Build_phylogeny.pl
#!/usr/bin/perl
use strict;
use warnings;

# phylogeny construction
# nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.fa.phy -n conca_sigpep -T 192 > tree.reports 2>&1 &

system("mkdir Phylogeny");
my @cmds;
my %hash1;
my $qua="Qualified_orth_recode.txt";
open QUA, $qua or die "can not open $qua\n";
while (<QUA>) {
        chomp;
        my @a=split;
        my ($Tnb, $Code);
        $Tnb =$a[2]+$a[3];
        $Code=$a[-1];
        $hash1{$a[-2]}++;

        if ($hash1{$a[-2]}==1) {
                my %hash2;
                my @b=split /\&\&/, $a[-2];
                foreach my $b (@b) {
                        $hash2{$b}++;
                }
                my $output="Orth_".$Tnb."spe_".$Code."code.phy";
                my $name  ="Orth_".$Tnb."spe_".$Code."code";
                open OUTPUT, ">Phylogeny/$output" or die "can not create $output\n"; 

                my $phy="conca_sigpep.fa.phy";
                open PHY, $phy or die "can not open $phy\n";
                while (<PHY>) {
                        chomp;
                        my @c=split;
                        if (/^\d+/) {
                                print OUTPUT "$Tnb  $c[1]\n";
                        } elsif ($hash2{$c[0]}) {
                                print OUTPUT "$_\n";
                        } else {
                                next;
                        }
                }
                my $cmd="raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s $output -n $name -T 190";
                push @cmds, $cmd;
        }
}

chdir("Phylogeny/");
foreach my $cmd (@cmds) {
        #print "$cmd\n";
        system($cmd);
}
```

```prepare_evo.pl
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

my $script_dir = dirname(abs_path($0));
my $qua="Qualified_orth_recode.txt";
open QUA, $qua or die "can not open $qua\n";
while (<QUA>) {
        chomp;
        my @a =split /\t/;
        my $id=$script_dir."/".$a[0]."_gene";
        my $dir="$id/OrthoFinder/Results_Dec14/Orthogroups";
        #print "$dir\n";
        chdir($dir);
        system("cp ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/prepare_input_paml.pl ./");
        my $cmd1="perl prepare_input_paml.pl --input orthologous_list_rep.txt ";
        $cmd1.="--seq_dir ~/Nocturnal_fish/Orthologous/pep/Phototransduction/sequences ";
        $cmd1.="--cor_list correlation.txt --output .";
        system($cmd1);
        chdir($script_dir);
}
```

```create_evo_tre.pl
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

my @phys=<RAxML_bestTree*>;
foreach my $phy (@phys) {
        (my $code)=$phy=~/RAxML_bestTree\.Orth_\d+spe_(\d+)code/;
        my $new="spe_".$code."code.tre";
        open NEW, ">$new\n" or die "can not create $new\n";
        open PHY, $phy or die "can not open $phy\n";
        while (<PHY>) {
                chomp;
                s/\:0\.[0-9]+//g;
                print NEW "$_\n";
        }
}
```

## Plot the sequence alignment
```bash
# OG0004341: CAHZ (PSGs && convergent evolution)
# kangjingliang@KangdeMacBook-Pro-2 二 12 30 2025 13:22:11 ~/Documents/2025/Nocturnal_fish/Target_genes/OG0004341
perl ../temp2.pl final_alignment_pep.fa > final_alignment_pep_sorted.fa

# OG0002364: GRK1 (PSGs && Relax selection)
# kangjingliang@KangdeMacBook-Pro-2 二 12 30 2025 13:22:11 ~/Documents/2025/Nocturnal_fish/Target_genes/OG0002364
perl ../temp2.pl final_alignment_pep.fa > final_alignment_pep_sorted.fa
```

## Functional enrichment analysis
### 1. Obtain the sequences of all orthologous genes for paml input 
```bash
# jlkang@hnu2024 Tue Sep 23 2025 10:53:25 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl create_orth_seq.pl > orth_seq.fasta
perl create_orth_paml.pl > final_orth_input_paml.fasta

# kangjingliang@KangdeMacBook-Pro-2 二 12 09 2025 14:47:03 ~/Documents/2025/Nocturnal_fish
cat Hyphy_PSGs/Hyphy_PSGs_genes.txt csubst/csubst_ccs_genes_id.txt Hyphy_relax/Intensified_genes.txt Hyphy_relax/Relaxed_genes.txt > Evolution_combine_genes.txt
```

## core clock genes
```bash
# kangjingliang@KangdeMacBook-Pro-2 二  9 23 2025 16:53:48 ~/Documents/2025/Nocturnal_fish
mkdir Core_clock_genes;cd Core_clock_genes
```

## Convergent site detection
```bash
# # ct discovery -a final_alignment_pep.fa -t ../config.tab -o discovery.output --fmt fasta
# h2076@h2076 Thu Sep 17 2026 11:14:45 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
# perl Run_caas.pl final_orth_input_paml.txt > run_caas.process 2>&1 &
# [1] 397584

# 只筛选pattern1的caas
perl Search_pattern1.pl > Total_ccas_pattern1.txt
```























