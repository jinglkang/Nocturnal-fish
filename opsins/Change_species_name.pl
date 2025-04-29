#!/usr/bin/perl
use strict;
use warnings;

my %hash=(
	'Zebrafish' => 'Danio_rerio',
	'Stickleback' => 'Gasterosteus_aculeatus',
	'Fugu' => 'Takifugu_rubripes',
	'Platyfish' => 'Xiphophorus_maculatus',
	'Medaka' => 'Oryzias_latipes',
	'Padel' => 'Pomacentrus_adelus',
	'Pmol' => 'Pomacentrus_moluccensis',
	'Apoly' => 'Acanthochromis_polyacanthus',
	'Acura' => 'Amblyglyphidodon_curacao',
	'Daru' => 'Dascyllus_aruanus',
	'Rgracilis' => 'Rhabdamia_gracilis',
	'Zviridiventer' => 'Zoramia_viridiventer',
	'Zleptacanthus' => 'Zoramia_leptacanthus',
	'Fthermalis' => 'Fibramia_thermalis',
	'Odoederleini' => 'Ostorhinchus_doederleini',
	'Ocookii' => 'Ostorhinchus_cookii',
	'Onovemfasciatus' => 'Ostorhinchus_novemfasciatus',
	'Onigrofasciatus' => 'Ostorhinchus_nigrofasciatus',
	'Onotatus' => 'Ostorhinchus_notatus',
	'Ocompressus' => 'Ostorhinchus_compressus',
	'Oangustatus' => 'Ostorhinchus_angustatus',
	'Ocyanosoma' => 'Ostorhinchus_cyanosoma',
	'Cquinquelineatus' => 'Cheilodipterus_quinquelineatus',
	'Cmacrodon' => 'Cheilodipterus_macrodon',
	'Cartus' => 'Cheilodipterus_artus',
	'Acrassiceps' => 'Apogon_crassiceps',
	'Fvariegata' => 'Fowleria_variegata',
	'Nfusca' => 'Nectamia_fusca',
	'Pmirifica' => 'Pterapogon_cf_mirifica',
	'Nviria' => 'Nectamia_viria',
	'Nsavayensis' => 'Nectamia_savayensis',
	'Amelas' => 'Apogonichthyoides_melas',
	'Abrevicaudatus' => 'Apogonichthyoides_brevicaudatus',
	'Snematoptera' => 'Sphaeramia_nematoptera',
	'Pexostigma' => 'Pristiapogon_exostigma',
	'Pfraenatus' => 'Pristiapogon_fraenatus',
	'Tfucata' => 'Taeniamia_fucata',
	'Tzosterophora' => 'Taeniamia_zosterophora'
	);

my $phy=$ARGV[0];
open PHY, $phy or die "can not open $phy\n";
while (<PHY>) {
	chomp;
	if (/^\d+/) {
		print "$_\n";
	} else {
		my @a=split;
		my $spe=$hash{$a[0]};
		print "$spe $a[1]\n";
	}
}
