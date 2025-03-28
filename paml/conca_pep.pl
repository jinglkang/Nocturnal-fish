#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
        chomp;
        my $orthdir=$_;
        chdir "$orthdir";
        my $cmd1="translateDna.pl -i final_alignment.fa > final_alignment_pep.fa";
        system($cmd1);
        chdir "/data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input";
}

my (%hash1, %hash2); my @spes;
my $orthpep="final_orth_input_paml.txt";
open ORTHPEP, $orthpep or die "can not open $orthpep\n";
while (<ORTHPEP>) {
        chomp;
        my $orthdir=$_;
        my $pep="$orthdir/final_alignment_pep.fa";
        my $spe;
        open PEP, $pep or die "can not open $pep\n";
        while (<PEP>) {
                chomp;
                if (/\>/) {
                        s/\>//;
                        $spe=$_;
                        $hash1{$spe}++;
                        push @spes, $spe if $hash1{$spe}==1;
                } else {
                        $hash2{$spe}.=$_;
                }
        }
}

foreach my $spe (@spes) {
        print ">$spe\n$hash2{$spe}\n";
}
