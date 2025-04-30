## False convergence and Random convergence detection
```bash
# use PER2 as an example, which have been detected with four convergent sites
# use codeml to estimate the branch lengths, amino acid frequencies and the best shape parameter for variable rates among sites (alpha) based on the amino acid sequences
# (base) jlkang@hnu2024 Wed Apr 30 2025 11:25:51 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0009073
vi spe.tre2
# ((Stickleback,((Fugu,((Platyfish,Medaka),(((Padel,Pmol),(Apoly,Acura)),Daru))),(((Rgracilis,((((Zviridiventer,Zleptacanthus),Fthermalis),((((Odoederleini,Ocookii),(Onovemfasciatus,Onigrofasciatus)),(Onotatus,((Ocompressus,(Oangustatus,Ocyanosoma)),Cquinquelineatus))),(Cmacrodon,Cartus))),(((Acrassiceps,Fvariegata),(((Nfusca,Pmirifica),(Nviria,Nsavayensis)),(Amelas,Abrevicaudatus))),(Snematoptera,(Pexostigma,Pfraenatus))))),Tfucata),Tzosterophora))),Zebrafish);
fasta2phy.pl final_alignment_pep.fa # phylip format: final_alignment_pep.fa.phy

```
