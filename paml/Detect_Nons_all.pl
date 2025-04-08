#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
        chomp;
        system("perl Detect_Nons.pl $_");
}
