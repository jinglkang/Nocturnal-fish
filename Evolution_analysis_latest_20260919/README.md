# Evolution of core clock genes (CCGs) in nocturnal fish
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
nohup perl prepare_input_paml_parallel.pl Target_orthologous_list_rep.txt >prepare_input_paml.process 2>&1 &
```

```Extract_target_orthogroups.pl
#!/usr/bin/perl
use strict;
use warnings;

my %dius=&build_dius();
my %gens=&build_genes();
my $anno="all_swissprot_diamond_ano_final.txt";
my %hash1;
open ANNO, $anno or die "can not open $anno\n";
while (<ANNO>) {
    chomp;
    my @a=split /\t/;
    (my $id)=$a[1]=~/sp\|.*\|(.*?)\_.*/;
    if (/^Zebrafish/ && $gens{$id}) {
        $hash1{$a[0]}=$id;
    }
}

my @headers;
my $orths="Orthogroups.tsv";
open ORTHS, $orths or die "can not open $orths\n";
while (<ORTHS>) {
    chomp;
    s/\,//g;
    my @a=split /\t/;
    if (/^Orthogroup/) {
        @headers=@a;
        print "GeneName\tDiuSpeNb\tNocSpeNb\tTotalSpeNb\t$_\n";
    } else {
        my ($nocspe, $diuspe);
        for (my $i = 1; $i < @a; $i++) {
            my @b=split /\s+/, $a[$i];
            if ($dius{$headers[$i]}) {
                $diuspe++ if @b>=1;
            } else {
                $nocspe++ if @b>=1;
            }
        }
        if ($diuspe && $nocspe) {
            my $Tspe=$diuspe+$nocspe;
            my @zebrs=split /\s+/, $a[-3];
            my $name1;
            if (@zebrs >= 1) {
                foreach my $id (@zebrs) {
                    if ($hash1{$id}) {
                        $name1=$hash1{$id};
                        last;
                    }
                }
            }
            if ($diuspe >=6 && $nocspe >= 10 && $Tspe <= 38 && $name1) {
                print "$name1\t$diuspe\t$nocspe\t$Tspe\t$_\n";
            }
        }
    }
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

sub build_genes {
    my %hash=(
        'OP1S1'=> 1,
        'OPSB'=> 1,
        'OPSD'=> 1,
        'OPSD1'=> 1,
        'OPSD2'=> 1,
        'OPSG'=> 1,
        'OPSG1'=> 1,
        'OPSG2'=> 1,
        'OPSG3'=> 1,
        'OPSG4'=> 1,
        'OPSR'=> 1,
        'OPSR1'=> 1,
        'OPSR2'=> 1,
        'OPSV'=> 1,
        'BHE40'=> 1,
        'BHE41'=> 1,
        'BMAL1'=> 1,
        'BMAL2'=> 1,
        'CIART'=> 1,
        'CLOCK'=> 1,
        'CRY1'=> 1,
        'CRY2'=> 1,
        'HLF'=> 1,
        'NFIL3'=> 1,
        'NPAS2'=> 1,
        'NR1D1'=> 1,
        'NR1D2'=> 1,
        'PER1'=> 1,
        'PER2'=> 1,
        'PER3'=> 1,
        'RORAA'=> 1,
        'RORAB'=> 1,
        'RORB'=> 1,
        'TEF'=> 1,
    );
    return(%hash);
}
```

```create_orth_rep.pl
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
my $subneed="Target_orthologous_list.txt";
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
        my $ID;
        ($orth{$sp}->{'ID'})?($ID=$orth{$sp}->{'ID'}):($ID="--");
        $info.=$ID."\t";
    }
    $info=~s/\s+$//;
    print "$info\n";
}
```









