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
        'Zebrafish'=> 1,
        'Stickleback'=> 1,
        'Fugu'=> 1,
        'Platyfish'=> 1,
        'Medaka'=> 1,
        'Daru'=> 1,
        'Apoly'=> 1,
        'Padel'=> 1,
        'Acura'=> 1,
        );

# compare the nonsynonymous position pep sequences one by one
my %hash1;
my @poss;
my @nocls=qw(Tzosterophora Ocompressus Ocyanosoma Oangustatus Onotatus Odoederleini Onovemfasciatus Onigrofasciatus Cartus Rgracilis Fthermalis Zviridiventer Zleptacanthus Snematoptera Pexostigma Pfraenatus Fvariegata Acrassiceps Nfusca Nsavayensis Nviria);
my @cleas=qw(Zebrafish Stickleback Fugu Platyfish Medaka Daru Apoly Padel Acura);
my @aspes=qw(Zebrafish Stickleback Fugu Platyfish Medaka Daru Apoly Padel Acura Tzosterophora Ocompressus Ocyanosoma Oangustatus Onotatus Odoederleini Onovemfasciatus Onigrofasciatus Cartus Rgracilis Fthermalis Zviridiventer Zleptacanthus Snematoptera Pexostigma Pfraenatus Fvariegata Acrassiceps Nfusca Nsavayensis Nviria);

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
                        # $hash2{$spepos}++;
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
