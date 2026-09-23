#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(getcwd);

my %dius=&build_dius();
my $orth="free_ratio_result.txt";
my %hash;
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
    chomp;
    my @a=split /\t/;
    if (/^Orth_id/) {
        next;
    } else {
        $dius{$a[1]}?($hash{$a[0]}->{'DIUS'}++):($hash{$a[0]}->{'NCOL'}++);
    }
}

my %hash1; my $header;
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
    chomp;
    my @a=split /\t/;
    my $file=$a[0]."_FreeRatio.txt";
    if (/^Orth_id/) {
        $header="$_\tType\n";
    } elsif ($hash{$a[0]}->{'DIUS'} && $hash{$a[0]}->{'NCOL'}) {
        if ($hash{$a[0]}->{'DIUS'} >= 6 && $hash{$a[0]}->{'NCOL'} >=10) {
            open FILE, ">>$file" or die "can not open $file\n";
            $hash1{$a[0]}++;
            print FILE "$header" if $hash1{$a[0]}==1;
            my $type;
            $dius{$a[1]}?($type="Diurnal"):($type="Nocturnal");
            print FILE "$_\t$type\n";
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
