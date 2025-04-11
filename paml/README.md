## Restart in HNU: nocturnal fish (Cardinalfishes)
### RUN Transdecoder
```bash
# download the refseq to annotate for opsin genes
# build the index for refseq
nohup wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/*protein.faa.gz &
zcat *gz >vertebrate_other.fa
# blast
makeblastdb -in vertebrate_other.fa -dbtype prot -input_type fasta
# diamond
diamond makedb --in vertebrate_other.fa -d vertebrate_other

# build the index for uniprot
# blast
makeblastdb -in uniprot-filtered-reviewed_yes.fasta -dbtype prot -input_type fasta -out uniprot
# diamond
diamond makedb --in uniprot-filtered-reviewed_yes.fasta -d uniprot
# jlkang@hnu2024 Mon Jan 13 16:11:30 /data2/jlkang/Nocturnal_fish
ll *.fasta|perl -alne '(my $nm)=$F[-1]=~/(.*)\.fasta/;$nm=~s/\_/ /g;print $nm' > species.txt
for i in *.fasta;do grep '>' ${i}|wc -l;done|less
```

```Run_transdecoder.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @cmds;
my @fasta=<*.fasta>;
foreach my $fas (@fasta) {
	(my $name)=$fas=~/(.*)\.fasta/;
	my $output=$name."_orf";
	my $cmd="TransDecoder.LongOrfs -t $fas -O $output";
	#print "$cmd\n";
	#system($cmd);
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(40);
foreach my $cmd (@cmds) {
	$manager->start and next;
	system($cmd);
	$manager->finish;
}
$manager -> wait_all_children;
```

```bash
# jlkang@hnu2024 Mon Jan 13 17:05:36 /data2/jlkang/Nocturnal_fish
vi Run_transdecoder.pl
nohup perl Run_transdecoder.pl >Run_transdecoder.process 2>&1 &
# [1] 12673
```
### Use the previous results
```bash
# 38 species: 28 nocturnal; 10 diurnal
# jlkang@hnu2024 Mon Jan 13 20:35:48 /data2/jlkang/Nocturnal_fish/Orthologous/pep
nohup orthofinder -f ./ > orthofinder.process 2>&1 &
# [1] 15314
# jlkang@hnu2024 Tue Jan 14 09:58:58 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan13/Orthogroups
less Orthogroups.GeneCount.tsv|head -n 2|perl -alne '$nb=@F;print $nb' # 40 coloumn; 38 species
# 3842 orthogroups with at least one sequence
less Orthogroups.GeneCount.tsv|perl -alne 'my ($i,$j);for ($i=1;$i<@F-1;$i++){$j++ if $F[$i]>=1};print "$F[0]" if $j==38' > orthologous_list.txt

# use "possvm" to divide the orthogroups in to smaller orthogroups
# install conda, possvm
# To activate this environment, use
#
#     $ conda activate ete3
#
# To deactivate an active environment, use
#
#     $ conda deactivate

# jlkang@hnu2024 Tue Jan 14 11:05:07 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan13
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ete3
conda install -c etetoolkit ete3 ete_toolchain
conda deactivate

# jlkang@hnu2024 Tue Jan 14 11:23:56 /data2/jlkang
mkdir software
# install conda in /data2/jlkang/software/
# create environment for possvm
conda create -n possvm
conda activate possvm
# install dependencies
conda install -c etetoolkit ete3==3.1.2
conda install -c bioconda pandas markov_clustering matplotlib numpy
conda install -c conda-forge networkx==3.0

# (possvm) jlkang@hnu2024 Thu Jan 16 15:12:14 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15
nohup perl build_sub_orth.pl > build_sub_orth.process 2>&1 &
# (base) jlkang@hnu2024 Thu Jan 16 19:05:13 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth
less sub_orth_genecount.txt|perl -alne 'my ($i,$j);for ($i=1;$i<@F;$i++){$j++ if $F[$i]>=1};print "$F[0]" if $j==38' > sub_ortho_38.txt
# 1326 # 1326 genes with at least one transcript in all of the 38 fish species
# only keep the longest transcript as the representative

# get the sequences belong to each suborth
# (base) jlkang@hnu2024 Thu Jan 23 10:16:27 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth
perl temp1.pl > sub_ortho_38seqs.txt
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

my @spes=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus 
    Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus 
	Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma 
	Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora 
	Zleptacanthus Zviridiventer Acura Apoly Daru Pmol Padel 
	Platyfish Fugu Medaka Stickleback Zebrafish);
my $head="suborth\t";
foreach my $spe (@spes) {
	$head.=$spe."\t";
}
$head=~s/\s+$//;
#print "$head\n";
my ($subinfo, $subfa);
my $subneed="sub_ortho_38.txt";
open SUBNEED, $subneed or die "can not open $subneed\n";
while (<SUBNEED>) {
	chomp;
	s/\s+$//;
	my $name=$_;
	my ($nm1, $nm2);
	my %orth;
	if (/\_/) {
		($nm1, $nm2)=$_=~/(.*)\_(.*)/;
		$subinfo="possvm_output/".$nm1."_tree.txt.ortholog_groups.csv";
		open SUBINFO, $subinfo or die "can not open $subinfo\n";
		while (<SUBINFO>) {
			chomp;
			my @a=split;
			if ($a[1] eq $nm2) {
				my ($spe, $geneid)=$a[0]=~/(.*?)\_(.*)/;
				$orth{$spe}.=$geneid.";";
			}		
		}			
	} else {
		$subfa="$name.fa";
		open SUBFA, $subfa or die "can not open $subfa\n";
		while (<SUBFA>) {
			chomp;
			if (/\>/) {
				s/\>//;
				(my $spe)=$_=~/(.*)\_/;
#				print "$_\t";
				$orth{$spe}=$_;
			}
		}
	}	
	my $info.=$name."\t";
	foreach my $spe (@spes) {
		my $info1=$orth{$spe};
		$info1=~s/\;$//;
		$info.=$info1."\t";
	}	
	$info=~s/\s+$//;
	print "$info\n";
}
```

```bash
# (base) jlkang@hnu2024 Thu Jan 23 10:19:26 /data2/jlkang/Nocturnal_fish/Orthologous/pep
cat *.fas > /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/all.fas
# only keep the longest transcript as the representative
# (base) jlkang@hnu2024 Thu Jan 23 10:54:39 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth
perl temp2.pl > sub_ortho_38seqs_rep.txt
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

my $allfas="all.fas";
my %len; my $name;
open ALLFAS, $allfas or die "can not open $allfas\n";
while (<ALLFAS>) {
        chomp;
        if (/^>/) {
                s/>//; $name=$_;
        } else {
                $len{$name}.=$_;
        }
}

my $list="sub_ortho_38seqs.txt";
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        if (/^suborth/) {
                print "$_\n";
        } else {
                my @a=split;
                my $info=$a[0]."\t";
                for (my $i = 1; $i < @a; $i++) {
                        my @b=split /\;/, $a[$i];
                        my %orth;
                        foreach $b (@b) {
                                my $seq=$len{$b};
                                if ($orth{$a[$i]}) {
                                        my ($old, $new);
                                        $old=$orth{$a[$i]}->{'len'};
                                        $new=length($len{$b});
                                        if ($new > $old) {
                                                $orth{$a[$i]}={
                                                        'len'=>length($len{$b}),
                                                        'ele'=>$b
                                                };
                                        }
                                } else {
                                        $orth{$a[$i]}={
                                                'len'=>length($len{$b}),
                                                'ele'=>$b
                                        };
                                }
                        }
                        $info.=$orth{$a[$i]}->{'ele'}."\t";
                }
                $info=~s/\s+$//;
                print "$info\n";
        }
}
```

```bash
# use the zebrafish sequences to annotate
# get the zebrafish Ensembl id
# (base) jlkang@hnu2024 Thu Jan 23 12:47:41 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth
less sub_ortho_38seqs_rep.txt|perl -alne 'next if /^sub/;my ($nm)=$F[-1]=~/.*\_(.*)/;print $nm' > sub_ortho_38seqs_rep_zebraID.txt
# kangjingliang@KangdeMacBook-Pro 四  1 23 12:49:44 ~/Documents/2025/Nocturnal_fish
scp jlkang@10.33.247.14:/data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/sub_ortho_38seqs_rep_zebraID.txt ./

# prepare the input for PAML
# (base) jlkang@hnu2024 Thu Jan 23 15:44:48 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth
mkdir paml_input; cp sub_ortho_38seqs_rep.txt paml_input/
# (base) jlkang@hnu2024 Thu Jan 23 16:20:22 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
less sub_ortho_38seqs_rep.txt|perl -alne 'my $info="$F[0]\t$F[-1]\t";for (my $i = 1; $i < @F-1; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;print $info' > ortho_list.txt
# (base) jlkang@hnu2024 Thu Jan 23 16:38:48 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
perl temp1.pl
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

my @peps=</data2/jlkang/Nocturnal_fish/Orthologous/pep/*.fas>;
my @nucs=</data2/jlkang/Nocturnal_fish/Orthologous/nuc/*.fas>;
foreach my $pep (@peps) {
        my $file=basename($pep);
        my ($name)=$file=~/(.*)\.fas/;
        my $new=$name."_pep.fasta";
        my $cmd="cp $pep $new";
        #print "$cmd\n";
        system($cmd);
}

foreach my $nuc (@nucs) {
        my $file=basename($nuc);
        my ($name)=$file=~/(.*)\.fas/;
        my $new=$name."_nuc.fasta";
        my $cmd="cp $nuc $new";
        #print "$cmd\n";
        system($cmd);
}
```

```bash
# (base) jlkang@hnu2024 Thu Jan 23 16:45:13 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
head -n 1 ortho_list.txt|perl -alne 'for (my $i = 1; $i < @F; $i++){print "$F[$i]\t$F[$i]_pep.fasta\t$F[$i]_nuc.fasta"}' > correlation.txt
# (base) jlkang@hnu2024 Thu Jan 23 16:47:41 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
nohup perl prepare_input_paml.pl --input ortho_list.txt --seq_dir . --cor_list correlation.txt --output . > prepare_input_paml.process 2>&1 &
# [1] 8262
# (base) jlkang@hnu2024 Sun Jan 26 12:28:40 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
wc -l final_orth_input_paml.txt # only 177 ortholgous genes can be used
# only 177 

# phylogenetic tree construction on these single orthlogous genes
# concatenate the protein sequences of these orthlogous genes
# 1. translate the nucleotide sequences to protein sequences (https://github.com/david-a-parry/translateDna/blob/master/translateDna.pl)
# translateDna.pl -i final_alignment.fa
# (base) jlkang@hnu2024 Mon Mar 17 20:46:09 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
perl conca_pep.pl >conca_sigpep.fa
# 40232 bp
```

```conca_pep.pl
#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
	chomp;
	my $orthdir=$_;
	chdir "$orthdir";
	my $cmd1="translateDna.pl -i final_alignment.fa > final_alignment_pep.fa";
	system($cmd1);
	chdir "/data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input";
}

my (%hash1, %hash2); my @spes;
my $orthpep="final_orth_input_paml.txt";
open ORTHPEP, $orthpep or die "can not open $orthpep\n";
while (<ORTHPEP>) {
	chomp;
	my $orthdir=$_;
	my $pep="$orthdir/final_alignment_pep.fa";
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
# construct the phylogenetic tree
# use the PROTGAMMAAUTO parameter to select the optimal amino acid substitution model, specify spotted gar as the outgroup, and evaluate the robustness of the result using 1000 bootstraps
# raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.phy -o Spottedgar -n conca_sigpep -T 192
# (base) jlkang@hnu2024 Mon Mar 17 23:34:18 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
fasta2phy.pl conca_sigpep.fa # generate phylip format: conca_sigpep.fa.phy
# we have 192 threads
# nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.fa.phy -o Spottedgar -n conca_sigpep -T 192 >tree.reports 2>&1 &
# (base) jlkang@hnu2024 Mon Mar 17 23:40:36 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/paml_input
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.fa.phy -n conca_sigpep -T 192 >tree.reports 2>&1 &
# [1] 297493
```

## Too few orhthologous genes, only select the longest transcripts as the represented genes for each species
```bash
# use all diurnal fishes as outgroup
# (base) jlkang@hnu2024 Fri Apr 11 21:47:37 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.fa.phy -o Zebrafish,Stickleback,Fugu,Platyfish,Medaka,Padel,Pmol,Apoly,Acura,Daru -n conca_sigpep_diurnal_outgroup -T 192 >tree.reports 2>&1 &
# [1] 198162
# use Zebrafish as outgroup, construct a new phylogeny
# (base) jlkang@hnu2024 Fri Apr 11 12:53:16 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input 
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.fa.phy -o Zebrafish -n conca_sigpep_zebrafish_root -T 192 >tree.reports 2>&1 &
# [1] 186740
# (base) jlkang@hnu2024 Tue Mar 18 15:55:08 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl temp1.pl > orthologous_list_rep.txt
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

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
my $subneed="orthologous_list.txt";
open SUBNEED, $subneed or die "can not open $subneed\n";
while (<SUBNEED>) {
	chomp;
	s/\s+$//; my $name=$_;
	my $fasta="../Orthogroup_Sequences/$name.fa";
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
	my $info="$name\t";
	foreach my $sp (@spes) {
		my $ID=$orth{$sp}->{'ID'};
		$info.=$ID."\t";
	}
	$info=~s/\s+$//;
	print "$info\n";
}
```

```bash
# (base) jlkang@hnu2024 Tue Mar 18 16:02:38 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
cp ../../sub_orth/paml_input/*.fasta ./
cp ../../sub_orth/paml_input/correlation.txt ./
cp orthologous_list_rep.txt paml_input/
# split orthologous_list_rep.txt for small files to run more quickly
# (base) jlkang@hnu2024 Tue Mar 18 16:29:28 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
split -l 40 orthologous_list_rep.txt smallfile
# (base) jlkang@hnu2024 Tue Mar 18 16:33:50 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup perl temp1.pl > prepare_input_paml.process 2>&1 &
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @cmds;
my @fasta=<smallfile*>;
foreach my $fas (@fasta) {
	my $cmd="perl prepare_input_paml.pl --input $fas --seq_dir . --cor_list correlation.txt --output .";
	#print "$cmd\n";
	#system($cmd);
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(100);
foreach my $cmd (@cmds) {
	$manager->start and next;
	system($cmd);
	$manager->finish;
}
$manager -> wait_all_children;
```

```bash
# (base) jlkang@hnu2024 Tue Mar 18 16:57:17 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup perl temp1.pl > prepare_input_paml.process 2>&1 & # make sure that the head order of "correlation.txt" is consistent with "orthologous_list_rep.txt"
# [1] 304019
# (base) jlkang@hnu2024 Wed Mar 19 10:00:20 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
wc -l final_orth_input_paml.txt # 1798 orthologous genes for the positive selection analysis
# concanated the sequences to build phylogenetic tree
# (base) jlkang@hnu2024 Wed Mar 19 10:20:49 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
perl conca_pep.pl >conca_sigpep.fa
# create the phylip format file conca_sigpep.fa.phy
fasta2phy.pl conca_sigpep.fa
# 30 species, 381363 amino acids 
# (base) jlkang@hnu2024 Wed Mar 19 10:20:49 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 1000 -s conca_sigpep.fa.phy -n conca_sigpep -T 192 > tree.reports 2>&1 &
# annotate all sequences
# swiss_pro_info.txt
# ref_pro_info.txt
# bioDBNet_db2db.pl
# (base) jlkang@hnu2024 Wed Mar 19 10:54:24 /data2/jlkang/software/database
perl extract_swiss_pro_info.pl >swiss_pro_info.txt
perl extract_refseq_pro_info.pl >ref_pro_info.txt
# (base) jlkang@hnu2024 Wed Mar 19 11:08:40 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
cp /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/sub_orth/all.fas ./
nohup annotate --fasta all.fas >annotate.reports 2>&1 & # run it after building the phylogenetic tree
# [2] 1058432

# estimate the evolutionary rate
# PAML: free-ratio
# (base) jlkang@hnu2024 Thu Mar 20 08:24:59 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
vi spe.tre
# (((((Snematoptera,(Pexostigma,Pfraenatus)),(((Tfucata,(Tzosterophora,((Zebrafish,Stickleback),(Fugu,((Platyfish,Medaka),(((Padel,Pmol),(Apoly,Acura)),Daru)))))),Rgracilis),(((Zviridiventer,Zleptacanthus),Fthermalis),((((Odoederleini,Ocookii),(Onovemfasciatus,Onigrofasciatus)),(Onotatus,((Ocompressus,(Oangustatus,Ocyanosoma)),Cquinquelineatus))),(Cmacrodon,Cartus))))),(Acrassiceps,Fvariegata)),((Nfusca,Pmirifica),(Nviria,Nsavayensis))),Amelas,Abrevicaudatus);
# (base) jlkang@hnu2024 Thu Mar 20 08:32:42 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup perl codeml_freeRatio_parallel.pl final_orth_input_paml.txt >free_ratio.process 2>&1 &
# [1] 1486788
# (base) jlkang@hnu2024 Fri Mar 28 22:48:54 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
perl extract_free_ratio.pl final_orth_input_paml.txt
```
### Free-ratio plot
```extract_free_ratio.pl
use strict;
use warnings;
open "fil", "$ARGV[0]";
my @names=<fil>;
my $header="Orth_id\tspecies\tlenth\t";
$header.="branch\tt\tN\tS\tdN/dS\tdN\tdS\tN*dN\tS*dS\n";
open BRANCH, ">>free_ratio_result.txt";
print BRANCH "$header";
foreach my $name (@names) {
        chomp $name;
        my $i;
        my ($len,$tag);
        my (%spe);
        open "fil1", "$name/free-ratio-result.txt";
        while (<fil1>) {
                chomp;
                $i++;
                my @a=split;
                if ($i==1) {
                        $len=$a[-1];
                } elsif (/^#/) {
                        $a[0]=~/#(\d+):/;
                        my $tag=$1;
                        $spe{$tag}=$a[-1];
                }elsif (@a==9 && $a[-1]=~/\d+\.\d+/) {
                        my $first=shift @a;
                        $first=~/\d+\.\.(\d+)/;
                        my $tag=$1;
                        $a[0]=~s/\.\./\_/;
                        my $len2=$a[2]+$a[3];
                        unless ($a[-3]>1 || $a[2]>$len || $len2>$len+50 || $a[-2]<1 || $a[-1]<1) {
                                print BRANCH "$name\t$spe{$tag}\t$len\t$_\n" if (exists $spe{$tag});
                        }
                }
        }
}
```

```bash
# kangjingliang@KangdeMacBook-Pro 五  3 28 22:53:09 ~/Documents/2025/Nocturnal_fish
perl extract_free_ratio.pl final_orth_input_paml.txt
scp jlkang@10.33.247.14:/data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/free_ratio_result.txt ./
# R plot
# ~/Documents/2025/Nocturnal_fish/free-ratio-plot.R
```
## Detect the convergence evolution in nocturnal fishes
```bash
# (base) jlkang@hnu2024 Tue Apr 08 20:24:21 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0013142
less final_alignment_pep.fa|grep '>'|perl -alne 's/\>//g;$info.=$_." ";END{print $info}'
# Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish
# (base) jlkang@hnu2024 Tue Apr 08 20:44:16 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
perl Detect_Nons_all.pl > convergent_evo_genes.txt
# kangjingliang@KangdeMacBook-Pro-2 二  4 08 21:08:55 ~/Documents/2025/Nocturnal_fish
scp all_swissprot_diamond_ano_final.txt jlkang@10.33.247.14:/data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input

# Orthogroups annotation
# (base) jlkang@hnu2024 Tue Apr 08 21:27:23 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
perl anno_orth.pl > convergent_evo_genes_ano.txt
```

```anno_orth.pl
#!/usr/bin/perl
use strict;
use warnings;

my $ann="all_swissprot_diamond_ano_final.txt";
my (%anno, %id);
open ANN, $ann or die "can not open $ann\n";
while (<ANN>) {
        chomp;
        my @a=split /\t/;
        if (/^Zebrafish/i) {
                $anno{$a[0]}=$a[1]."\t".$a[-1];
        }
}

my $list="orthologous_list_rep.txt";
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        next if /^Orth/i;
        my @a=split /\t/;
        $id{$a[0]}=$a[-1];
}

my $cove="paml_input/convergent_evo_genes.txt";
open COVE, $cove or die "can not open $cove\n";
while (<COVE>) {
        chomp;
        my @a=split /\t/;
        my $zeb=$id{$a[0]};
        my $an =$anno{$zeb};
        print "$a[0]\t$an\t$a[1]\t$a[2]\n";
}
```

## Detect the positive selected genes in the ancestoral node of all nocturnal fish species
```bash
# (base) jlkang@hnu2024 Wed Apr 09 13:41:09 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
vi spe_nocAces.tre
# (Zebrafish,(Stickleback,((Fugu,((Medaka,Platyfish),(Daru,((Acura,Apoly),(Pmol,Padel))))),(Tzosterophora,(Tfucata,(Rgracilis,(((Snematoptera,(Pfraenatus,Pexostigma)),((Fvariegata,Acrassiceps),((Amelas,Abrevicaudatus),((Nsavayensis,Nviria),(Pmirifica,Nfusca))))),((Fthermalis,(Zleptacanthus,Zviridiventer)),((Cartus,Cmacrodon),(((Onigrofasciatus,Onovemfasciatus),(Ocookii,Odoederleini)),(Onotatus,(Cquinquelineatus,(Ocompressus,(Ocyanosoma,Oangustatus)))))))))))#1)));
# perl codeml.pl --input temp/$temp --model branch-site --dir . --output_suf Nocturnal --tree spe_nocAces.tre --icode 0 --omega 1.2
vi codeml_parallel_PSGs.pl
```

```codeml_parallel_PSGs.pl
use strict;
use warnings;
use Parallel::ForkManager;

# prepare_input_paml.pl
my $fm=$ARGV[0];
system("mkdir temp");

my @cmds;
open FM, $fm or die "can not open $fm\n";
while (<FM>) {
        chomp; my @a=split;
        my $temp=$a[0].".txt";
        open TEMP, ">temp/$temp" or die "can not create >temp/$temp\n";
        print TEMP "$_\n";
        my $cmd="perl codeml.pl --input temp/$temp --model branch-site --dir . --output_suf Nocturnal --tree spe_nocAces.tre --icode 0 --omega 1.2";
        push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(100);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;
system("rm -rf temp/");
```

```bash
# (base) jlkang@hnu2024 Wed Apr 09 13:54:57 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup perl codeml_parallel_PSGs.pl final_orth_input_paml.txt >codeml.process 2>&1 &
# [1] 154795
```
