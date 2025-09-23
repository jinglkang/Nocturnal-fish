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
                $anno{$a[0]}=$a[1];
        }
}

# ~/software/database/swiss_pro_info.txt
my $uni="/home/jlkang/software/database/uniprot_sprot.fasta";
my %seq; my ($ID);
open UNI, $uni or die "can not open $uni\n";
while (<UNI>) {
        chomp;
        if (/>/) {
                s/\>//;
                my @a=split;
                $ID=$a[0];
        } else {
                $seq{$ID}.=$_;
        }
}

my $list="orthologous_list_rep.txt";
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        next if /^Orth/i;
        my @a=split /\t/;
        if ($anno{$a[-1]}) {
                my $uni=$anno{$a[-1]};
                print ">$a[0]\n$seq{$uni}\n" if $seq{$uni};
        }
}
