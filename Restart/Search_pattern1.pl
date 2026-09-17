#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

my $header="Gene\tTrait\tPosition\tSubstitution\tPvalue\tPattern\tFFGN\tFBGN\tGFG\tGBG\tMFG\tMBG\tFFG\tFBG\tMS";
print "$header\n";
my @caas=<*_caas.txt>;
foreach my $caa (@caas) {
        (my $orth)=$caa=~/(.*)_caas\.txt/;
        open CAA, $caa or die "can not open $caa\n";
        while (<CAA>) {
                chomp;
                my @a=split /\t/;
                if (/pattern1/) {
                        my $info;
                        for (my $i = 1; $i < @a; $i++) {
                                $info.=$a[$i]."\t";
                        }
                        $info=~s/\s+$//;
                        print "$orth\t$info\n";
                }
        }
}
