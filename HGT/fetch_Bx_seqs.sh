module load Anaconda3/2022.05

#### Files and directories
wd=/mnt/parscratch/users/bo1lpg/NEOF-project/bx-cluster/bx-new-assembly
query_list=/mnt/parscratch/users/bo1lpg/NEOF-project/bx-cluster/evolutionary-history-bx-cluster/approach1-miniprot/target_genes.txt
query_location=/mnt/parscratch/users/bo1lpg/NEOF-project/bx-cluster/trees-bx-genes/query
DB_location=/mnt/parscratch/users/bo1lpg/Genome_dbs
DB_to_use=${wd}/species.txt
blast_results=${wd}/blast_results
consensus=/mnt/parscratch/users/bo1lpg/NEOF-project/phylobased-LGT/nested_scripts/consensus.pl

## Build blast databases for selected species
cd ${wd}
source activate blast
mkdir databases
while read genome; do
  makeblastdb -in ${DB_location}/${genome}/${genome}_helixer.cds -dbtype nucl -out databases/${genome}_cds
done < ${DB_to_use}

## Run blastn with relaxed parameters for the maize genes and filter by match length = 300
mkdir -p relaxed_parameters/blast_results
while read gene; do
  while read genome; do
    blastn -query ${query_location}/${gene}.fa -db databases/${genome}_cds \
      -word_size 9 -gapopen 3 -penalty -2 -gapextend 2 \
      -outfmt '6 sseqid length sseq pident evalue' | awk '$2 >= 300' > relaxed_parameters/blast_results/${gene}_${genome}_filtered.txt
    done < ${DB_to_use}
done < ${query_list}

#### Select matches (only hit sequence) and put them in a fasta file --> one per gene
cd ${wd}
mkdir files_to_align
while read gene; do
  while read genome; do
    cat blast_results/${gene}_${genome}_filtered.txt | while read match; do
      echo ${match} | awk -v species=${genome} '{print ">" species "." $1 "\n" $3}' >> files_to_align/${gene}.fa;
    done
  done < ${DB_to_use}
  sed -i 's/-//g' files_to_align/${gene}.fa
done < ${query_list}

#### Save the number of matches per gene
ls files_to_align | while read file ; do grep ">" files_to_align/"$file" \
    | sort | uniq | wc -l | sed 's/^/'$file'\t/g' >> blast_matches_per_gene.txt ; done

#### Align the blast fragments to the original sequence using MAFFT & unwrap alingment
cd ${query_location}
source activate mafft
mkdir ${wd}/ALN1-mafft
unset MAFFT_BINARIES

while read line ; do
  mafft --thread 8 --localpair --maxiterate 1000 --addfragments ${wd}/files_to_align/"$line".fa ${query_location}/${line}.fa > ${wd}/ALN1-mafft/"$line"
done < ${query_list}

while read line ; do
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' \
    < ${wd}/ALN1-mafft/"$line" > ${wd}/ALN1-mafft/"$line"_unwrap ;
done < ${query_list}

while read line ; do
  tail -n +2 ${wd}/ALN1-mafft/"$line"_unwrap > ${wd}/ALN1-mafft/"$line" ; rm ${wd}/ALN1-mafft/"$line"_unwrap
done < ${query_list}
source deactivate

#### Identify blast matches represented by only one sequence fragment and move these to their own folder
mkdir ${wd}/ALN2-unique
while read line ; do
  grep ">" ${wd}/ALN1-mafft/"$line" | cut -f 2 -d ">" | sort | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | grep "\s1$" | cut -f 1 -d ' ' | while read line2 ; do
    grep "$line2$" -A 1 ${wd}/ALN1-mafft/"$line" >> ${wd}/ALN2-unique/"$line"
  done
done < ${query_list}

#### Identify blast matches represented by more than one sequence fragment and move these to their own folder
mkdir ${wd}/ALN3-duplicated
while read line ; do
  grep ">" ${wd}/ALN1-mafft/"$line" | cut -f 2 -d ">" | sort | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | grep -v "\s1$" | cut -f 1 -d ' ' | while read line2 ; do
    grep "$line2$" -A 1 ${wd}/ALN1-mafft/"$line" >> ${wd}/ALN3-duplicated/"$line"_"$line2"
  done
done < ${query_list}

#### Generate a single consensus sequence for blast matches represented by more than one sequence
mkdir ${wd}/ALN4-consensus
module load BioPerl/1.7.8-GCCcore-12.2.0
ls ${wd}/ALN3-duplicated | while read line ; do perl ${consensus} -in ${wd}/ALN3-duplicated/"$line"  -out ${wd}/ALN4-consensus/"$line" -iupac; done
ls ${wd}/ALN4-consensus  | while read line ; do sed -i '/>/c\>'$line'' ${wd}/ALN4-consensus/"$line" ; done
while read line ; do
  sed -i 's/>'$line'_/>/g' ${wd}/ALN4-consensus/"$line"_*
done < ${query_list}

#### Merge the consensus and unique sequences into a single alignment
mkdir ${wd}/Fasta_mafft_alignments
while read line ; do
  cat ${wd}/ALN2-unique/"$line" ${wd}/ALN4-consensus/"$line"_* > ${wd}/Fasta_mafft_alignments/"$line"
done < ${query_list}

while read line ; do
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${wd}/Fasta_mafft_alignments/"$line" > ${wd}/Fasta_mafft_alignments/"$line"_unwrap
done < ${query_list}

while read line ; do
  tail -n +2 ${wd}/Fasta_mafft_alignments/"$line"_unwrap > ${wd}/Fasta_mafft_alignments/"$line" ; rm ${wd}/Fasta_mafft_alignments/"$line"_unwrap
done < ${query_list}

#### Check whether the number of sequences in the alignment is correct - the 'check' files must be manually inspected
ls find -name '${wd}/Fasta_mafft_alignments/*' -size 0 -delete
ls ${wd}/Fasta_mafft_alignments | while read line ; do grep ">" ${wd}/Fasta_mafft_alignments/"$line" | sort | uniq | wc -l | sed 's/^/'$line'\t/g' >> ${wd}/sequences_in_alignment.txt ; done
paste ${wd}/sequences_in_alignment.txt ${wd}/blast_matches_per_gene.txt | awk 'BEGIN { OFS = "\t" } NR == 0 { $5 = "diff." } NR >= 0 { $5 = $2 - ($4+1) } 1' \
    | cut -f 1,5 | grep -v "\s0$" > ${wd}/Check_for_errors.txt
cat ${wd}/sequences_in_alignment.txt | grep -v "\s0" > ${wd}/sequences_in_alignment_no0.txt
cat ${wd}/sequences_in_alignment_no0.txt ${wd}/blast_matches_per_gene.txt | cut -f 1 | sort | uniq -u > ${wd}/missing_alignments.txt
cp ${wd}/blast_matches_per_gene.txt ${wd}/blast_matches_per_gene_no-missing-aln.txt
cat ${wd}/missing_alignments.txt | while read line ; do sed -i '/'$line'/d' ${wd}/blast_matches_per_gene_no-missing-aln.txt ; done
paste ${wd}/sequences_in_alignment_no0.txt ${wd}/blast_matches_per_gene_no-missing-aln.txt \
    | awk 'BEGIN { OFS = "\t" } NR == 0 { $5 = "diff." } NR >= 0 { $5 = $2 - ($4+1) } 1' | cut -f 1,5 | grep -v "\s0$" > ${wd}/Check_for_errors2.txt

#### Trim alignments with trimAl
cd ${wd}
mkdir trimmed_aln
source activate trimal
while read line; do
  trimal -in Fasta_mafft_alignments/${line} -out trimmed_aln/${line}_auto -automated1
done < ${query_list}
