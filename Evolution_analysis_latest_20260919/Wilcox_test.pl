#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(getcwd);

my @ratios=<*FreeRatio.txt>;
foreach my $ratio (@ratios) {
        (my $orth)=$ratio=~/(.*)\_FreeRatio\.txt/;
        my $rscrip="Rscript.R";
        open RSCRIP, ">$rscrip" or die "can not create $rscrip\n";
        my $cmd1="Mdata <- read.table(\"$ratio\",header = TRUE)";
        my $cmd2="Noct<-Mdata[Mdata[,\"Type\"]==\"Nocturnal\",]";
        my $cmd3="Diur<-Mdata[Mdata[,\"Type\"]==\"Diurnal\",]";
        my $cmd4="wilcox.test(Noct\$dN.dS,Diur\$dN.dS,exact = FALSE)";
        print RSCRIP "$cmd1\n$cmd2\n$cmd3\n$cmd4\n";
        my $output=`Rscript $rscrip`;
        my $result=$orth."_wilcox_result.txt";
        open RESULT, ">$result" or die "can not create $result\n";
        print RESULT "$output";
}
