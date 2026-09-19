# Evolution of core clock genes (CCGs) and opsins in nocturnal fish
## Extract the orthogroups if any of zebrafish genes was annotated as target genes
### Target genes: opsins, core clock genes
```bash
# h2076@h2076 Fri Sep 18 2026 22:58:26 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups
# 要求orthogroup里面必须有斑马鱼的序列，且该序列被注释为目标基因的名字，同时夜行鱼数量大于6，日行鱼数量大于10
perl Extract_target_orthogroups.pl > Target_Orthogroups.txt
# 包含所有38个物种序列的orthogroups已经被检查过了，所以此次关注小于38个物种的orthogroups
less Target_Orthogroups.txt|perl -alne 'next if /Orthogroup/i;my $nb=$F[3];print $F[4] if $nb<38' > Target_orthologous_list.txt
perl create_orth_rep.pl > Target_orthologous_list_rep.txt
mkdir Target_genes; cd Target_genes
# h2076@h2076 Fri Sep 18 2026 23:33:09 ~/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/Target_genes
cp ../paml_input/*.fasta ./; cp ../paml_input/correlation.txt ./; cp ../Target_orthologous_list_rep.txt ./
cp ../paml_input/prepare_input_paml.pl ./
# perl prepare_input_paml.pl --input Target_orthologous_list_rep.txt --seq_dir . --cor_list correlation.txt --output .
nohup perl prepare_input_paml_parallel.pl Target_orthologous_list_rep.txt >prepare_input_paml.process 2>&1 &
```
