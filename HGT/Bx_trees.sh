module load Anaconda3/2022.05
source activate iqtree

i=$(expr $SLURM_ARRAY_TASK_ID)

#### Directories and input files
wd=/mnt/parscratch/users/bo1lpg/NEOF-project/bx-cluster/bx-new-assembly/default_parameters/trimmed_aln
fasta_list=/mnt/parscratch/users/bo1lpg/NEOF-project/bx-cluster/evolutionary-history-bx-cluster/approach1-miniprot/target_genes.txt
outdir=${wd}/trees

#### Scripts
fasta_to_phylip=/mnt/parscratch/users/bo1lpg/Aristidoideae/trees/Fasta2Phylip.pl

#### Step 2: create directories and convert fasta to phylip
head -$i ${fasta_list} | tail -1 | while read line ; do mkdir -p ${outdir}/individual/"$line" ; mkdir -p ${outdir}/combined ; mkdir -p ${outdir}/logs ;  \
  cd ${outdir}/individual/"$line" ; perl ${fasta_to_phylip} ${wd}/"$line"_auto "$line" ; done

#### Step 3: run sms and make a copy of the tree
head -$i ${fasta_list} | tail -1 | while read line ; do cd ${outdir}/individual/"$line" ; iqtree -s "$line" -m MFP -B 1000 -bnni -alrt 1000 -T AUTO -ntmax 14; cp *.contree ${outdir}/combined ; done

#### Step 4: clean up by moving log files
mv slurm-*_$i.out ${outdir}/logs
