#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

# Detect_Nons.pl
my %seq; my $name;
my $orth=$ARGV[0];
my $fas="$orth/final_alignment_pep.fa";
open FAS, $fas or die "can not open $fas\n";
while (<FAS>) {
        chomp;
        if (/^>/) {
                s/\>//;
                $name=$_;
        } else {
                $seq{$name}.=$_;
        }
}

my %cleaner=(
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

# compare the nonsynonymous position pep sequences one by one
my %hash1;
my @poss;
my @nocls=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer);
my @cleas=qw(Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish);
my @aspes=qw(Abrevicaudatus Acrassiceps Amelas Cartus Cmacrodon Cquinquelineatus Fthermalis Fvariegata Nfusca Nsavayensis Nviria Oangustatus Ocompressus Ocookii Ocyanosoma Odoederleini Onigrofasciatus Onotatus Onovemfasciatus Pexostigma Pfraenatus Pmirifica Rgracilis Snematoptera Tfucata Tzosterophora Zleptacanthus Zviridiventer Acura Apoly Daru Pmol Padel Platyfish Fugu Medaka Stickleback Zebrafish);

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
                        $hash2{$spepos}++;
                        push @cleas_pos, $spepos;
                }
                foreach my $spe (@nocls) {
                        my $spepos=$hash1{$spe}->{$pos};
                        push @nocls_pos, $spepos;
                }
                my @isect = intersect(@cleas_pos, @nocls_pos);
                my $numb=keys %hash2;
                unless (@isect) {
                        print "$orth\t$numb\t$info\n";
                }
        }
}
