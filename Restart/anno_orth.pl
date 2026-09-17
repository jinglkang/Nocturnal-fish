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

my $cove="Total_ccas_pattern1.txt";
open COVE, $cove or die "can not open $cove\n";
while (<COVE>) {
        chomp;
        my @a=split /\t/;
        print "UniprotID\tGene_description\t$_\n" if /^Gene/;
        my $zeb=$id{$a[0]};
        my $an =$anno{$zeb};
        print "$an\t$_\n";
}
