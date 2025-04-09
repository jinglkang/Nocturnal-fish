#!/usr/bin/perl
use strict;
use warnings;

my %ANO;
open ANO, $ARGV[0] or die "can not open $ARGV[0]\n";
while (<ANO>) {
        chomp;
        my @a=split /\t/;
        $ANO{$a[0]}=$a[1]."\t".$a[2];
}

my @psgs=<postively_selected_genes/*.txt>;
foreach my $psg (@psgs) {
        my ($spe, $ge)=$psg=~/postively_selected_genes\/(.*)-psg-(.*)\.txt/;
        open PSG, $psg or die "can not open $psg\n";
        while (<PSG>) {
                chomp;
                my @a=split /\t/;
                my $info;
                for (my $i = 1; $i < @a; $i++) {
                        $info.=$a[$i]."\t";
                }
                $info=~s/\s+$//;
                print "$spe\t$ge\t$ANO{$ge}\t$info\n";
        }
}
