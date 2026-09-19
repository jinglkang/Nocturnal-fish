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
