module load Anaconda3/2022.05

#### Directories and files
wd=/mnt/parscratch/users/bo1lpg/NEOF-project/sequencing/2025-assembly
species='Zuloagaea_bulbosa'
hifi=${wd}/raw/hifi/hifi_reads.fastq
assembly=${wd}/A06-Helixer_annotation/Zuloagaea_bulbosa_assembly.fa
organelle_refs=${wd}/A09_MarkOrganelles/ref_organelles

cd ${wd}/A09_MarkOrganelles

#### Map HiFi reads to assembly with minimap
source activate minimap2
minimap2 -t 16 -ax map-hifi ${assembly} ${hifi} | \
  samtools sort -o aln.bam
samtools index aln.bam

#### Compute per-contig coverage
samtools coverage aln.bam | \
  awk 'NR>1 {print $1, $7}' OFS="\t" \
  > contig_coverage.txt
source deactivate

#### Blast known chloroplast genomes against assembly
source activate blast
makeblastdb -in ${assembly} -dbtype nucl -out assembly

blastn \
  -query ${organelle_refs}/chloroplasts.fasta \
  -db assembly \
  -out cp_vs_genome.tsv \
  -evalue 1e-20 -perc_identity 90 -num_threads 16 \
  -outfmt "6 qseqid sseqid length pident qstart qend sstart send bitscore"

blastn \
  -query ${organelle_refs}/mitochondria.fasta \
  -db assembly \
  -out mt_vs_genome.tsv \
  -evalue 1e-20 -perc_identity 90 -num_threads 16 \
  -outfmt "6 qseqid sseqid length pident qstart qend sstart send bitscore"
source deactivate

#### Collapse blast hits per contig
source activate bedtools
awk '{
  s = ($7 < $8 ? $7 : $8);
  e = ($7 < $8 ? $8 : $7);
  print $2, s-1, e
}' OFS="\t" cp_vs_genome.tsv > cp_hits.bed

awk '{
  s = ($7 < $8 ? $7 : $8);
  e = ($7 < $8 ? $8 : $7);
  print $2, s-1, e
}' OFS="\t" mt_vs_genome.tsv > mt_hits.bed

bedtools sort -i cp_hits.bed | \
bedtools merge -i - | \
awk '{sum[$1] += ($3 - $2)} END {for (c in sum) print c, sum[c]}' \
  > cp_aligned_bp_per_contig.txt

bedtools sort -i mt_hits.bed | \
bedtools merge -i - | \
awk '{sum[$1] += ($3 - $2)} END {for (c in sum) print c, sum[c]}' \
  > mt_aligned_bp_per_contig.txt

#### Get contig length
cut -f1,2 ${assembly}.fai > contig_lengths.txt

#### Ratio between collapsed blast hit and contig length
awk '
NR==FNR {len[$1]=$2; next}
{
  if ($1 in len)
    print $1, $2, len[$1], $2/len[$1]
}
' contig_lengths.txt cp_aligned_bp_per_contig.txt \
  > cp_fraction_per_contig.txt

awk '
NR==FNR {len[$1]=$2; next}
{
  if ($1 in len)
    print $1, $2, len[$1], $2/len[$1]
}
' contig_lengths.txt mt_aligned_bp_per_contig.txt \
> mt_fraction_per_contig.txt


#### Join with coverage data
awk '
NR==FNR {cov[$1]=$2; next}
{
  print $1, $2, $3, $4, cov[$1]
}
' contig_coverage.txt cp_fraction_per_contig.txt \
> cp_fraction_with_coverage.txt

awk '
NR==FNR {cov[$1]=$2; next}
{
  print $1, $2, $3, $4, cov[$1]
}
' contig_coverage.txt mt_fraction_per_contig.txt \
> mt_fraction_with_coverage.txt
