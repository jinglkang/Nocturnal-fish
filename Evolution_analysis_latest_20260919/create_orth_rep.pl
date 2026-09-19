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
