wd=/mnt/parscratch/users/bo1lpg/NEOF-project/LGT_new_assembly/
cds=${wd}/CDS/Zuloagaea_bulbosa_01/Zuloagaea_bulbosa_01_final.cds.fa
aln_folder=${wd}/results_08_mark_dups/Fasta_mafft_alignments/selected-aln
gff=/mnt/parscratch/users/bo1lpg/NEOF-project/sequencing/2025-assembly/A06-Helixer_annotation/Zuloagaea_bulbosa.gff3

## Make Directory
cd ${wd}
mkdir LGT_function_blast
cd LGT_function_blast

## Get query list
ls ${aln_folder} | sed "s/__/:/g" | cut -d ":" -f 2 > query.list


while read gene; do
  grep ${gene} ${gff} >> annotations_alldups.txt
done < dup.list

## Get query sequences
while read gene; do
  grep -A 1 ${gene} ${cds} >> sequences.fa
done < gene_IDs_Jan.txt

## Make database with swissprot
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
source activate blast
makeblastdb -in uniprot_sprot.fasta -dbtype prot

## Blast CDS against swissprot
blastx -query sequences.fa -db uniprot_sprot.fasta -max_target_seqs 1 \
  -outfmt "6 qseqid sallseqid sallgi qstart qend sstart send evalue bitscore salltitles" > LGT_candidates_blast_results.txt

blastx -query sequences.fa -db ${wd}/LGT_function_blast/uniprot_sprot.fasta -max_target_seqs 1 \
  -outfmt "6 qseqid sallseqid sallgi qstart qend sstart send evalue bitscore salltitles" > other_genes_blast_results.txt

## Keep only one blast hit per sequence (higher bitscore)
awk '
{
    q = $1
    bs = $9

    if (!(q in seen)) {
        order[++n] = q
        seen[q] = 1
    }

    if (!(q in best_bs) || bs > best_bs[q]) {
        best_bs[q] = bs
        best_line[q] = $0
    }
}
END {
    for (i = 1; i <= n; i++) print best_line[order[i]]
}
' LGT_candidates_blast_results.txt > LGT_candidates_blast_results_filtered.txt
