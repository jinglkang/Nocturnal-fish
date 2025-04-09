#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

# Detect_Nons.pl
my %seq; my $name;
my $fas="final_alignment_pep.fa";
open FAS, $fas or die "can not open $fas\n";
while (<FAS>) {
        chomp;
        if (/^>/) {
                s/\>//g;
                s/\_.*//g;
                $name=$_;
        } else {
                $seq{$name}.=$_;
        }
}

my %cleaner=(
        'Medaka'=> 1,
        'Platyfish'=> 1,
        'Daru'=> 1,
        'Padel'=> 1,
        'Pmol'=> 1,
        );

# compare the nonsynonymous position pep sequences one by one
my %hash1;
my @poss;
my @nocls=qw(Tfucata Rgracilis Pmirifica Pfraenatus Pexostigma Fvariegata Acrassiceps Amelas Abrevicaudatus Nfusca Nsavayensis Nviria Cmacrodon Cartus Onigrofasciatus Onovemfasciatus Ocookii Odoederleini Onotatus Ocompressus Oangustatus Ocyanosoma);
my @cleas=qw(Medaka Platyfish Daru Padel Pmol);
my @aspes=qw(Medaka Platyfish Daru Padel Pmol Tfucata Rgracilis Pmirifica Pfraenatus Pexostigma Fvariegata Acrassiceps Amelas Abrevicaudatus Nfusca Nsavayensis Nviria Cmacrodon Cartus Onigrofasciatus Onovemfasciatus Ocookii Odoederleini Onotatus Ocompressus Oangustatus Ocyanosoma);

&Build_pos_hash(\@aspes);

sub Build_pos_hash {
        my ($grp)=@_;
        my @grp=@{$grp};
        my $len;
        foreach my $spe (@grp) {
                my $seq=$seq{$spe};
                $len=length($seq);
                for (my $i = 0; $i < $len; $i++) {
                        my $spepos=substr($seq,$i,1);
                        $hash1{$spe}->{$i}=$spepos;
                }
        }

        for (my $i = 0; $i < $len; $i++) {
                my %hash2;
                my $pos=$i;
                my $newp=$pos+1;
                my $info=$newp.":";
                foreach my $spe (@aspes) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $info.=$spe."($spepos);";
                }

                my (@cleas_pos, @nocls_pos);
                foreach my $spe (@cleas) {
                        my $spepos=$hash1{$spe}->{$pos};
                        push @cleas_pos, $spepos;
                }
                foreach my $spe (@nocls) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $hash2{$spepos}++;
                        push @nocls_pos, $spepos;
                }
                my @isect = intersect(@cleas_pos, @nocls_pos);
                my $numb=keys %hash2;
                unless (@isect) {
                        print "OPSB\t$numb\t$info\n";
                }
        }
}
