# Nocturnal-fish
## Select the fasta files needed
```bash
# 28 nocturnal fish species
# Kang@fishlab3 Tue Oct 17 15:07:12 ~/下载/Cardinalfish_assembled_transcripts
cp M06_P_mirifica_assembled_transcripts.fasta Needed_fasta/Pterapogon_cf_mirifica.fasta
cp M07_N_fusca_assembled_transcripts.fasta Needed_fasta/Nectamia_fusca.fasta
cp M122_Nsavayensis_assembled_transcripts.fasta Needed_fasta/Nectamia_savayensis.fasta
cp M83_N_vir_assembled_transcripts.fasta Needed_fasta/Nectamia_viria.fasta
cp M191_Abrevicaudatus_assembled_transcripts.fasta Needed_fasta/Apogonichthyoides_brevicaudatus.fasta
cp M10_Amelas_assembled_transcripts.fasta Needed_fasta/Apogonichthyoides_melas.fasta
cp M97_S_nem_assembled_transcripts.fasta Needed_fasta/Sphaeramia_nematoptera.fasta
cp M156_Acrassiceps_assembled_transcripts.fasta Needed_fasta/Apogon_crassiceps.fasta
cp M163_Fvariegata_assembled_transcripts.fasta Needed_fasta/Fowleria_variegata.fasta
cp M127_Pfraenatus_assembled_transcripts.fasta Needed_fasta/Pristiapogon_fraenatus.fasta
cp M127_Pfraenatus_assembled_transcripts.fasta Needed_fasta/Pristiapogon_exostigma.fasta
cp M04_Rgracilis_assembled_transcripts.fasta Needed_fasta/Rhabdamia_gracilis.fasta
cp M135_Fthermalis_assembled_transcripts.fasta Needed_fasta/Fibramia_thermalis.fasta
cp E8A_Zviridiventer_assembled_transcripts.fasta Needed_fasta/Zoramia_viridiventer.fasta
cp M17_Zleptacanthus_assembled_transcripts.fasta Needed_fasta/Zoramia_leptacanthus.fasta
cp M128_Cartus_assembled_transcripts.fasta Needed_fasta/Cheilodipterus_artus.fasta
cp M001_Cquinquelineatus_assembled_transcripts.fasta Needed_fasta/Cheilodipterus_quinquelineatus.fasta
cp M141_Cmacrodon_assembled_transcripts.fasta Needed_fasta/Cheilodipterus_macrodon.fasta
cp M148_Tzosterophora_assembled_transcripts.fasta Needed_fasta/Taeniamia_zosterophora.fasta
cp M03_T_fucata_assembled_transcripts.fasta Needed_fasta/Taeniamia_fucata.fasta
cp M004_O_cyanosoma_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_cyanosoma.fasta
cp M003_Oangustatus_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_angustatus.fasta
cp M123_Onigrofasciatus_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_nigrofasciatus.fasta
cp M119_Onovemfasciatus_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_novemfasciatus.fasta
cp M05_O_cookii_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_cookii.fasta
cp M002_Ocompressus_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_compressus.fasta
cp M02_O_doederleini_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_doederleini.fasta
cp M12_O_notatus_assembled_transcripts.fasta Needed_fasta/Ostorhinchus_notatus.fasta
```
## Prepare the input files for orthlogous genes detection
```bash
# Kang@fishlab3 Tue Oct 17 15:47:35 /media/HDD
mkdir Nocturnal_fish_project
# Kang@fishlab3 Tue Oct 17 16:40:12 /media/HDD/Nocturnal_fish_project
cp ~/下载/Cardinalfish_assembled_transcripts/Needed_fasta/*.fasta ./
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my @fasta=<*.fasta>;
foreach my $fas (@fasta) {
	(my $name)=$fas=~/(.*)\.fasta/;
	my $output=$name."_orf";
	my $cmd="TransDecoder.LongOrfs -t $fas -O $output";
	#print "$cmd\n";
	system($cmd);
}
```

```temp3.pl
#!/usr/bin/perl
use strict;
use warnings;

my $nucdir="Orthologous_genes/nuc";
my $pepdir="Orthologous_genes/pep";
system("mkdir -p $nucdir");
system("mkdir -p $pepdir");
my @dirs=<*_orf>;
foreach my $dir (@dirs) {
        my ($nm1, $nm2)=$dir=~/(\D{1}).*\_(.*?)\_orf/;
        my $name=$nm1.$nm2;
        my $nuc =$dir."/longest_orfs.cds";
        my $pep =$dir."/longest_orfs.pep";
        my $nucnew="$nucdir/$name".".fas";
        my $pepnew="$pepdir/$name".".fas";
        my $i;
        open NUC, $nuc or die "can not open $nuc\n";
        open NUCN,">$nucnew" or die "can not create $nucnew\n";
        while (<NUC>) {
                chomp;
                if (/>/) {
                        $i++;
                        my $nameH=$name."_".$i;
                        print NUCN ">$nameH\n";
                } else {
                        print NUCN "$_\n";
                }
        }

        $i=0;
        open PEP, $pep or die "can not open $pep\n";
        open PEPN,">$pepnew" or die "can not create $pepnew\n";
        while (<PEP>) {
                chomp;
                if (/>/) {
                        $i++;
                        my $nameH=$name."_".$i;
                        print PEPN ">$nameH\n";
                } else {
                        print PEPN "$_\n";
                }
        }
}
```

```bash
# detect ORF
# Kang@fishlab3 Thu Oct 19 10:33:48 /media/HDD/Nocturnal_fish_project
nohup perl temp1.pl > detect_orf.process 2>&1 &
# [1] 18267
# Kang@fishlab3 Thu Oct 19 15:08:29 /media/HDD/Nocturnal_fish_project
perl temp3.pl
# Kang@fishlab3 Thu Oct 19 14:48:14 /media/HDD/white_island/Compevo/orth16_new
cp Zebrafish.cds.fasta Stickleback.cds.fasta Platyfish.cds.fasta Medaka.cds.fasta Fugu.cds.fasta Apoly.cds.fasta Daru.cds.fasta Padel.cds.fasta Pmol.cds.fasta Acura.cds.fasta /media/HDD/Nocturnal_fish_project/Orthologous_genes/nuc
cp Zebrafish.pep.fasta Stickleback.pep.fasta Platyfish.pep.fasta Medaka.pep.fasta Fugu.pep.fasta Apoly.pep.fasta Daru.pep.fasta Padel.pep.fasta Pmol.pep.fasta Acura.pep.fasta /media/HDD/Nocturnal_fish_project/Orthologous_genes/pep
```

```/media/HDD/Nocturnal_fish_project/Orthologous_genes/pep/temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my @fass=<*.pep.fasta>;
foreach my $fas (@fass) {
	(my $name)=$fas=~/(.*)\.pep\.fasta/;
	my $name1=$name.".fas";
	system("mv $fas $name1");
}
```

```bash
# Kang@fishlab3 Thu Oct 19 14:58:19 /media/HDD/Nocturnal_fish_project/Orthologous_genes/pep
perl temp1.pl
```

```/media/HDD/Nocturnal_fish_project/Orthologous_genes/nuc/temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my @fass=<*.cds.fasta>;
foreach my $fas (@fass) {
	(my $name)=$fas=~/(.*)\.cds\.fasta/;
	my $name1=$name.".fas";
	system("mv $fas $name1");
}
```

```bash
# Kang@fishlab3 Thu Oct 19 15:02:06 /media/HDD/Nocturnal_fish_project/Orthologous_genes/nuc
perl temp1.pl
```

```bash
# orthologous genes detection
# Kang@fishlab3 Thu Oct 19 15:14:12 /media/HDD/Nocturnal_fish_project/Orthologous_genes/pep
nohup orthofinder -f ./ >orthofinder.process 2>&1 &
# [1] 2234

# swissprot to annotate all genes
# Kang@fishlab3 Fri Oct 20 16:31:38 /media/HDD/Nocturnal_fish_project/Orthologous_genes
cat pep/*.fas > all.fas
# Transfer to HPC to annotate
# jlkang@hpc2021 Fri Oct 20 16:34:46 /lustre1/g/sbs_schunter/Kang
mkdir Nocturnal_fish
# Kang@fishlab3 Fri Oct 20 16:35:49 /media/HDD/Nocturnal_fish_project/Orthologous_genes
scp all.fas jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Nocturnal_fish
# Transfer the database data
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 20 16:43:16 ~
scp -r swiss-prot jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang

# diamond blastp -q all.fasta -e 1e-5 --sensitive -k 1 -d /lustre1/g/sbs_schunter/Kang/swiss-prot/uniprot-filtered-reviewed_yes.fasta --out all_swissprot_diamond_ano.txt  > diamond_blastp.process
module load diamond
```

```script_diamond.cmd
#!/bin/bash
#SBATCH --job-name=diamond        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
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
diamond blastp -q all.fasta -e 1e-5 --sensitive -k 1 -d /lustre1/g/sbs_schunter/Kang/swiss-prot/uniprot-filtered-reviewed_yes.fasta --out all_swissprot_diamond_ano.txt  > diamond_blastp.process
# print the end time
```

```bash
# jlkang@hpc2021 Fri Oct 20 16:47:20 /lustre1/g/sbs_schunter/Kang/Nocturnal_fish
sbatch script_diamond.cmd
```
