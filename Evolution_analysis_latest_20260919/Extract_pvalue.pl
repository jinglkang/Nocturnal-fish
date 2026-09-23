#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(getcwd);

my %dius=&build_dius();
my @ratios=<*FreeRatio.txt>;
print "OrthId\tMnoc\tRelat\tMdiu\tPvalue\n";
foreach my $ratio (@ratios) {
        (my $orth)=$ratio=~/(.*)\_FreeRatio\.txt/;

        (my $wilc)=$orth."_wilcox_result.txt";
        my $pval;
        open WILC, "$wilc" or die "can not open $wilc\n";
        while (<WILC>) {
                chomp;
                if (/p-value\s+=\s+(.*)/) {
                        $pval=$1;
                }
        }

        my ($noc, $diu, $i, $j);
        open ORTH, $ratio or die "can not open $ratio\n";
        while (<ORTH>) {
                chomp;
                my @a=split;
                next if /^Orth_id/;
                if ($dius{$a[1]}) {
                        $i++;
                        $diu+=$a[7];
                } else {
                        $j++;
                        $noc+=$a[7];
                }
        }
        my ($Mdiu,$Mnoc);
        $Mdiu=$diu/$i;
        $Mnoc=$noc/$j;

        if ($Mnoc > $Mdiu) {
                print "$orth\t$Mnoc\t>\t$Mdiu\t$pval\n";
        } elsif ($Mnoc < $Mdiu) {
                print "$orth\t$Mnoc\t<\t$Mdiu\t$pval\n";
        } else {
                print "$orth\t$Mnoc\t=\t$Mdiu\t$pval\n";
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
