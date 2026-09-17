#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

my $list=$ARGV[0]; # The list
my $outd="CAAS";
unless (-d $outd) {
    mkdir $outd;
}

# ct discovery -a final_alignment_pep.fa -t ../config.tab -o discovery.output --fmt fasta

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
    chomp;
    my @a=split;
    my $orth=$a[0];
    my $alig="$orth/final_alignment_pep.fa"; # the alignment
    my $outl="$outd/$orth"."_caas.txt";
    my $cmd ="ct discovery -a $alig -t config.tab -o $outl --fmt fasta";
#   print "$cmd\n";
    push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(80);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
