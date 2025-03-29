### Grep opsin genes for each species
```bash
# kangjingliang@KangdeMacBook-Pro 三  3 26 09:12:42 ~/Documents/2025/Nocturnal_fish/Gene_filter
less Zebrafish_gene_description.txt|grep 'opsin'|grep -v 'teleost\|retinal\|exorh\|vertebrate\|parietopsin\|parapinopsin\|tissue\|like\|mRNA'|grep -i 'cone\|RHO'|perl -alne '@a=split /\t/;print if $a[1]' > Zebrafish_opsins.txt
less Fugu_gene_description.txt|grep 'opsin'|grep -v 'teleost\|retinal\|exorh\|vertebrate\|parietopsin\|parapinopsin\|tissue\|like\|mRNA'|grep -i 'cone\|RHO'|perl -alne '@a=split /\t/;print if $a[1]' > Fugu_opsins.txt
less Medaka_gene_description.txt|grep 'opsin'|grep -v 'teleost\|retinal\|exorh\|vertebrate\|parietopsin\|parapinopsin\|tissue\|like\|mRNA'|grep -i 'cone\|RHO'|perl -alne '@a=split /\t/;print if $a[1]' > Medaka_opsins.txt
less Stickleback_gene_description.txt|grep 'opsin'|grep -v 'teleost\|retinal\|exorh\|vertebrate\|parietopsin\|parapinopsin\|tissue\|like\|mRNA'|grep -i 'cone\|RHO'|perl -alne '@a=split /\t/;print if $a[1]' > Stickleback_opsins.txt
less Platyfish_gene_description.txt|grep 'opsin'|grep -v 'teleost\|retinal\|exorh\|vertebrate\|parietopsin\|parapinopsin\|tissue\|like\|mRNA'|grep -i 'cone\|RHO'|perl -alne '@a=split /\t/;print if $a[1]' > Platyfish_opsins.txt
```
```bash
# HKU HPC
# romeodan@hpc2021 Fri Oct 20 16:47:20 /lustre1/g/sbs_schunter/Kang/Nocturnal_fish
sbatch script_diamond.cmd
# Submitted batch job 1975816
```

```script_diamond.cmd
#!/bin/bash
#SBATCH --job-name=diamond        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jinglkang@outlook.com     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --time=7-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err

# print the start time
date
diamond blastp -q all.fas -e 1e-5 --sensitive -k 1 -d /lustre1/g/sbs_schunter/Kang/swiss-prot/uniprot-filtered-reviewed_yes.fasta --out all_swissprot_diamond_ano.txt  > diamond_blastp.process
# print the end time
```

### HPC environment is not correct after i change "./bashrc"
```bash
#PATH="$HOME/.local/bin:$HOME/bin:$PATH"
export PATH=/home/romeodan/.local/bin:/home/romeodan/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/:/lustre1/g/sbs_schunter/Kang/software/ncbi-blast-2.16.0+/bin:/lustre1/g/sbs_schunter/Kang/software/bin
PATH="/lustre1/g/sbs_schunter/Kang/software/ncbi-blast-2.16.0+/bin:/lustre1/g/sbs_schunter/Kang/software/bin:$PATH"
```

### Annotate and detect the target genes
```bash
# (base) kang1234@celia-PowerEdge-T640 Wed Mar 26 19:58:33 ~/swiss-prot
nohup diamond blastp -q all.fas -e 1e-5 --sensitive -k 1 -d uniprot-filtered-reviewed_yes.fasta --out all_swissprot_diamond_ano.txt > diamond_blastp.process 2>&1 &
# [1] 5413
# (base) kang1234@celia-PowerEdge-T640 Wed Mar 26 20:39:30 ~/swiss-prot
perl temp1.pl > all_swissprot_diamond_ano_final.txt

# kangjingliang@KangdeMacBook-Pro 三  3 26 23:29:58 ~/Documents/2025/Nocturnal_fish
#########################
# get the putative Opsins
mkdir Opsins
less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OPSR/i' > Opsins/OPSR.txt
less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OPSV/i' > Opsins/OPSV.txt
less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OP1S1/i' > Opsins/OP1S1.txt
less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OPSB/i' > Opsins/OPSB.txt
less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OPSD/i' > Opsins/OPSD.txt
less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OPSG/i' > Opsins/OPSG.txt

less all_swissprot_diamond_ano_final.txt|perl -alne 'my @a=split /\|/, $F[1];print if $a[-1]=~/OPSR/i || $a[-1]=~/OPSV/i || $a[-1]=~/OP1S1/i || $a[-1]=~/OPSB/i || $a[-1]=~/OPSD/i || $a[-1]=~/OPSG/i' > Opsins_blast.txt 
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my $bla=$ARGV[0];
my %hash;
open BLA, $bla or die "can not open $bla\n";
while (<BLA>) {
        chomp;
        my @a=split /\t/;
        my ($spe, $id, $len);
        $id=$a[0]; $len=$a[4];
        ($spe)=$a[0]=~/(.*?)\_/;
        if ($hash{$spe}) {
                my $oldlen=$hash{$spe}->{'LEN'};
                my $newlen=$len;
                if ($newlen > $oldlen) {
                        $hash{$spe}={
                        'ID'  => $id,
                        'LEN' => $len,
                        'INF' => $_
                        };
                }
        } else {
                $hash{$spe}={
                        'ID'  => $id,
                        'LEN' => $len,
                        'INF' => $_
                        };
        }
}

foreach my $spe (sort keys %hash) {
        my $info=$hash{$spe}->{'INF'};
        print "$info\n";
}
```

```bash
# kangjingliang@KangdeMacBook-Pro 四  3 27 00:38:02 ~/Documents/2025/Nocturnal_fish/Opsins
perl temp1.pl OPSB.txt > OPSB_target.txt
perl temp1.pl OPSD.txt > OPSD_target.txt
perl temp1.pl OPSG.txt > OPSG_target.txt
perl temp1.pl OPSR.txt > OPSR_target.txt

# kangjingliang@KangdeMacBook-Pro 四  3 27 00:44:28 ~/Documents/2025/Nocturnal_fish/Opsins
for i in *.txt;do wc -l ${i};done
#       1 OP1S1.txt
#      77 OPSB.txt
#      34 OPSB_target.txt
#     125 OPSD.txt
#      37 OPSD_target.txt
#      98 OPSG.txt
#      38 OPSG_target.txt
#      41 OPSR.txt
#      25 OPSR_target.txt
#       5 OPSV.txt

less OPSB_target.txt|perl -alne '$hash{$F[3]}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}' # select align length >=346
# 103	1
# 318	1
# 346	1
# 352	31
less OPSB_target.txt|perl -alne 'print $F[0] if $F[3] >=346' > OPSB_target_id.txt

less OPSD_target.txt|perl -alne '$hash{$F[3]}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}' # select align length >=329
# 114	2
# 116	2
# 137	1
# 144	1
# 188	1
# 227	1
# 228	1
# 289	1
# 329	1
# 330	11
# 331	5
# 353	8
# 354	2
less OPSD_target.txt|perl -alne 'print $F[0] if $F[3] >=329' > OPSD_target_id.txt

less OPSG_target.txt|perl -alne '$hash{$F[3]}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}' # select align length >=334
# 114	1
# 179	1
# 185	2
# 186	2
# 187	1
# 296	1
# 334	2
# 343	27
# 349	1
less OPSG_target.txt|perl -alne 'print $F[0] if $F[3] >=334' > OPSG_target_id.txt

less OPSR_target.txt|perl -alne '$hash{$F[3]}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}' # select align length >=356
# 189	1
# 356	1
# 357	23
less OPSR_target.txt|perl -alne 'print $F[0] if $F[3] >=356' > OPSR_target_id.txt

# kangjingliang@KangdeMacBook-Pro 四  3 27 11:14:29 ~/Documents/2025/Nocturnal_fish/Opsins
scp *_target_id.txt kang1234@10.64.139.91:~/nocturnal_fish
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 27 11:16:44 ~/nocturnal_fish
mkdir Opsins; mv *_target_id.txt Opsins/
perl temp1.pl *_target_id.txt # get the target sequences to fasta file: OPSB_target_id.fas; OPSD_target_id.fas; OPSG_target_id.fas; OPSR_target_id.fas
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%hash1, %hash2); my ($id1, $id2);
my $pep="../all_pep.fas";
open PEP, $pep or die "can not open $pep\n";
while (<PEP>) {
	chomp;
	if (/\>/) {
		s/\>//;
		$id1=$_;
	} else {
		$hash1{$id1}.=$_;
	}
}

my $nuc="../all_nuc.fas";
open NUC, $nuc or die "can not open $nuc\n";
while (<NUC>) {
	chomp;
	if (/\>/) {
		s/\>//;
		$id2=$_;
	} else {
		$hash2{$id2}.=$_;
	}
}

my @tar=@ARGV;
foreach my $tar (@tar) {
	(my $name)=$tar=~/(.*)\.txt/;
	my $output1=$name."_pep.fas";
	my $output2=$name."_nuc.fas";
	open OUTPUT1, ">$output1" or die "can not create $output1\n";
	open OUTPUT2, ">$output2" or die "can not create $output2\n";
	open TAR, $tar or die "can not open $tar\n";
	while (<TAR>) {
		chomp;
		my $seq1=$hash1{$_};
		my $seq2=$hash2{$_};
		print OUTPUT1 ">$_\n$seq1\n";
		print OUTPUT2 ">$_\n$seq2\n";
	}
	close TAR;
	close OUTPUT1;
	close OUTPUT2;
}
```

### align opsin genes to construct phylogenetic tree and codeml
```bash
# get the nucleotide first
# Kang@fishlab3 Thu Mar 27 11:48:08 /media/HDD/Nocturnal_fish_project/Orthologous_genes/nuc
cat *.fas > all_nuc.fas
scp all_nuc.fas kang1234@10.64.139.91:~/nocturnal_fish
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 27 12:52:22 ~/nocturnal_fish/Opsins
mkdir OPSB;mv OPSB_target_id_pep.fas OPSB_target_id_nuc.fas OPSB/
mkdir OPSD;mv OPSD_target_id_pep.fas OPSD_target_id_nuc.fas OPSD/
mkdir OPSG;mv OPSG_target_id_pep.fas OPSG_target_id_nuc.fas OPSG/
mkdir OPSR;mv OPSR_target_id_pep.fas OPSR_target_id_nuc.fas OPSR/
# for paml input and phylogenetic tree: align; remove gap; translate

# (base) kang1234@celia-PowerEdge-T640 Thu Mar 27 13:28:08 ~/nocturnal_fish/Opsins/OPSB 
clustalo -i OPSB_target_id_pep.fas -t Protein -o pep.aln --outfmt=fa
pal2nal.pl pep.aln OPSB_target_id_nuc.fas -output fasta >pal2nal.fasta
Gblocks pal2nal.fasta -b4=10 -b5=n -b3=5 -t=c
perl temp1.pl # output align nucleotide sequences
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
fasta2phy.pl final_alignment_pep.fa > final_alignment_pep.phy
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s final_alignment_pep.phy -n conca_sigpep -T 24 >tree.reports 2>&1 &
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
open FIL2, "pal2nal.fasta-gb" or die "can not open pal2nal.fasta-gb\n";
open FIL3, ">final_alignment.phy" or die "can not create final_alignment.phy\n";
open FIL4, ">final_alignment.fa" or die "can not create final_alignment.fa\n";
my $SPE; my (%hash1, %hash2); my @spe;
while (<FIL2>) {
    chomp;
    if (/>/) {
            s/>//; $SPE=$_; push @spe, $SPE;
    } else {
            s/ //g; s/\d/-/g;
            $hash1{$SPE}.=$_;
    }
}
my $i=0;
foreach my $spe (@spe) {
    my $seq=$hash1{$spe};
    for (my $i = 0; $i < length($seq); $i+=3) {
            my $codon=substr($seq,$i,3);
            if ($codon=~/taa/i||$codon=~/tga/i||$codon=~/tag/i) {
                    $codon="---"; # replace stop codon to gap
            }
            $hash2{$spe}.=$codon;
    }
    $i++; if ($i==1) {
            my $num=@spe; my $len=length($hash2{$spe});
            print FIL3 "$num  $len\n";
            print FIL3 "$spe  $hash2{$spe}\n";
            print FIL4 ">$spe\n$hash2{$spe}\n";
    } else {
            print FIL3 "$spe  $hash2{$spe}\n";
            print FIL4 ">$spe\n$hash2{$spe}\n";
    }
}
```

```bash
# OPCB
# kangjingliang@KangdeMacBook-Pro 四  3 27 13:35:34 ~/Documents/2025/Nocturnal_fish/Opsins
mkdir OPSB; cd OPSB
# kangjingliang@KangdeMacBook-Pro 四  3 27 13:36:16 ~/Documents/2025/Nocturnal_fish/Opsins/OPSB
scp kang1234@10.64.139.91:~/nocturnal_fish/Opsins/OPSB/RAxML* ./

# OPSD
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 27 13:39:55 ~/nocturnal_fish/Opsins/OPSD
cp ../OPSB/temp1.pl ./
clustalo -i OPSD_target_id_pep.fas -t Protein -o pep.aln --outfmt=fa
pal2nal.pl pep.aln OPSD_target_id_nuc.fas -output fasta >pal2nal.fasta
Gblocks pal2nal.fasta -b4=10 -b5=n -b3=5 -t=c
perl temp1.pl # output align nucleotide sequences
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
fasta2phy.pl final_alignment_pep.fa > final_alignment_pep.phy
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s final_alignment_pep.phy -n conca_sigpep -T 24 >tree.reports 2>&1 &
# [1] 19080
# kangjingliang@KangdeMacBook-Pro 四  3 27 13:35:34 ~/Documents/2025/Nocturnal_fish/Opsins
mkdir OPSD; cd OPSD
scp kang1234@10.64.139.91:~/nocturnal_fish/Opsins/OPSD/RAxML* ./

# OPSG
cp ../OPSB/temp1.pl ./
clustalo -i OPSG_target_id_pep.fas -t Protein -o pep.aln --outfmt=fa
pal2nal.pl pep.aln OPSG_target_id_nuc.fas -output fasta >pal2nal.fasta
Gblocks pal2nal.fasta -b4=10 -b5=n -b3=5 -t=c
perl temp1.pl # output align nucleotide sequences
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
fasta2phy.pl final_alignment_pep.fa > final_alignment_pep.phy
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s final_alignment_pep.phy -n conca_sigpep -T 24 >tree.reports 2>&1 &
# [1] 20569
# kangjingliang@KangdeMacBook-Pro 四  3 27 13:35:34 ~/Documents/2025/Nocturnal_fish/Opsins
mkdir OPSG; cd OPSG
scp kang1234@10.64.139.91:~/nocturnal_fish/Opsins/OPSG/RAxML* ./

# OPSR
cp ../OPSB/temp1.pl ./
clustalo -i OPSR_target_id_pep.fas -t Protein -o pep.aln --outfmt=fa
pal2nal.pl pep.aln OPSR_target_id_nuc.fas -output fasta >pal2nal.fasta
Gblocks pal2nal.fasta -b4=10 -b5=n -b3=5 -t=c
perl temp1.pl # output align nucleotide sequences
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
fasta2phy.pl final_alignment_pep.fa > final_alignment_pep.phy
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s final_alignment_pep.phy -n conca_sigpep -T 24 >tree.reports 2>&1 &
# [1] 21496
# kangjingliang@KangdeMacBook-Pro 四  3 27 13:35:34 ~/Documents/2025/Nocturnal_fish/Opsins
mkdir OPSR; cd OPSR
scp kang1234@10.64.139.91:~/nocturnal_fish/Opsins/OPSR/RAxML* ./

# select the orthologous genes of corresponsding species for each opsin to construct a phylogenetic tree
# concatenate the orthologous protein sequences of corresponsding species
# (base) jlkang@hnu2024 Sat Mar 29 08:02:08 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
perl conca_spePEP.pl OPSB_target_id.txt > OPSB_target_spe.fasta
fasta2phy.pl OPSB_target_spe.fasta # output OPSB_target_spe.phy
perl conca_spePEP.pl OPSB_target_id.txt > OPSB_target_spe.fasta
fasta2phy.pl OPSB_target_spe.fasta # output OPSB_target_spe.phy

perl conca_spePEP.pl OPSD_target_id.txt > OPSD_target_spe.fasta
fasta2phy.pl OPSD_target_spe.fasta # output OPSD_target_spe.phy

perl conca_spePEP.pl OPSG_target_id.txt > OPSG_target_spe.fasta
fasta2phy.pl OPSG_target_spe.fasta # output OPSG_target_spe.phy

perl conca_spePEP.pl OPSR_target_id.txt > OPSR_target_spe.fasta
fasta2phy.pl OPSR_target_spe.fasta # output OPSR_target_spe.phy

# construct the phylogenetci tree
# opsin_phy.sh
raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s OPSB_target_spe.phy -n OPSB -T 192
raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s OPSD_target_spe.phy -n OPSD -T 192
raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s OPSG_target_spe.phy -n OPSG -T 192
raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s OPSR_target_spe.phy -n OPSR -T 192
# (base) jlkang@hnu2024 Sat Mar 29 08:10:37 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input
nohup sh opsin_phy.sh > tree.reports 2>&1 &
# [1] 38653
```
