## False convergence and Random convergence detection
### 1. codeml estimation branch length etc.
```bash
# use PER2 as an example, which have been detected with four convergent sites
# 1. use codeml to estimate the branch lengths, amino acid frequencies and the best shape parameter for variable rates among sites (alpha) based on the amino acid sequences
# (base) jlkang@hnu2024 Wed Apr 30 2025 11:25:51 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0009073
vi spe.tre2
# ((Stickleback,((Fugu,((Platyfish,Medaka),(((Padel,Pmol),(Apoly,Acura)),Daru))),(((Rgracilis,((((Zviridiventer,Zleptacanthus),Fthermalis),((((Odoederleini,Ocookii),(Onovemfasciatus,Onigrofasciatus)),(Onotatus,((Ocompressus,(Oangustatus,Ocyanosoma)),Cquinquelineatus))),(Cmacrodon,Cartus))),(((Acrassiceps,Fvariegata),(((Nfusca,Pmirifica),(Nviria,Nsavayensis)),(Amelas,Abrevicaudatus))),(Snematoptera,(Pexostigma,Pfraenatus))))),Tfucata),Tzosterophora))),Zebrafish);
fasta2phy.pl final_alignment_pep.fa # phylip format: final_alignment_pep.fa.phy
codeml estimate_before_stimulation.ctr # it doesn't matter if occurs "error: end of tree file."
```

### 2. evolver for amino acid sequence simulation
#### Based on "mlc" (result file)
```branch_length
# branch_length
((Stickleback: 0.641521, ((Fugu: 0.131312, ((Platyfish: 0.102181, Medaka: 0.112472): 0.048615, (((Padel: 0.000004, Pmol: 0.004373): 0.008782, (Apoly: 0.026817, Acura:
 0.013352): 0.000004): 0.021468, Daru: 0.105587): 0.037622): 0.014017): 0.001577, (((Rgracilis: 0.026783, ((((Zviridiventer: 0.000004, Zleptacanthus: 0.000004): 0.008
970, Fthermalis: 0.009050): 0.013498, ((((Odoederleini: 0.004457, Ocookii: 0.008918): 0.000004, (Onovemfasciatus: 0.000004, Onigrofasciatus: 0.004487): 0.004476): 0.0
00004, (Onotatus: 0.009010, ((Ocompressus: 0.008968, (Oangustatus: 0.000004, Ocyanosoma: 0.013470): 0.000004): 0.004454, Cquinquelineatus: 0.004479): 0.000004): 0.000
004): 0.000004, (Cmacrodon: 0.000004, Cartus: 0.000004): 0.008970): 0.008894): 0.000004, (((Acrassiceps: 0.037325, Fvariegata: 0.027598): 0.003944, (((Nfusca: 0.00443
8, Pmirifica: 0.000004): 0.000004, (Nviria: 0.008901, Nsavayensis: 0.000004): 0.000004): 0.000004, (Amelas: 0.004441, Abrevicaudatus: 0.008882): 0.000004): 0.000004):
 0.000004, (Snematoptera: 0.031735, (Pexostigma: 0.000004, Pfraenatus: 0.000004): 0.004561): 0.004306): 0.000004): 0.000004): 0.017994, Tfucata: 0.004460): 0.000004,
Tzosterophora: 0.000004): 0.075914): 0.040696): 0.400904, Zebrafish: 0.400022);
```

```alpha_gamma
# alpha_gamma
.42248 3
```

```amino_acid_freq
# amino_acid_freq
0.05293 0.03954 0.03937 0.04472 0.02530 0.04915 0.06258 0.05187 0.02117 0.05284 0.07478 0.06279 0.02419 0.03758 0.06534 0.13332 0.05399
 0.00489 0.03405 0.06959
A R N D C Q E G H I L K M F P S T W Y V
```

```bash
# create "MCaa.dat" in the current directory and run "evolver" which will use MCaa automatically
# (base) jlkang@hnu2024 Wed Apr 30 2025 20:24:53 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0009073
evolver 7
ls -lt # list the files according the date (the lastest to earliest)
# the out put is: ancestral.txt; mc.txt; siterates.txt
# And then estimate the ratio of convergence and non-convergence in the selected site in the 1000 replicates
```
