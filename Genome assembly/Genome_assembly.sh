
module load Anaconda3/2022.05
source activate assembly

#### Directories and files
wd=/mnt/parscratch/users/bo1lpg/NEOF-project/sequencing/2025-assembly
hifi=${wd}/raw/hifi
hic=${wd}/raw/hic
species='Zuloagaea_bulbosa'

#### Generate fastq from bam file for hifi raw data
cd ${hifi}
bam2fastq multiple_movies.hifi_reads.bam -o hifi_reads.fastq

#### Run hifiasm
cd ${wd}
mkdir A01-PhasedAssembly
cd A01-PhasedAssembly
hifiasm -o ${species}_hifi_hi-c_n2.asm --n-hap 2 -t 48 \
  --h1 ${hic}/Zbulbosa-1A_1.fq.gz \
  --h2 ${hic}/Zbulbosa-1A_2.fq.gz \
  ${hifi}/hifi_reads.fastq

# Convert GFA to FASTA
cd ${wd}/A01-PhasedAssembly/
gfatools gfa2fa ${species}_hifi_hi-c_n2.asm.hic.p_utg.gfa > ${species}_hifi_hi-c_n2.asm.hic.p_utg.fa

# Basic stats
seqkit stats ${species}_hifi_hi-c_n2.asm.hic.p_utg.fa > ../unitig_stats.txt

# Filter by length (keep >=10kb; change as needed)
seqkit seq -m 10000 ${species}_hifi_hi-c_n2.asm.hic.p_utg.fa -o ${species}_hifi_hi-c_n2.asm.hic.p_utg.10kb.fa
seqkit stats ${species}_hifi_hi-c_n2.asm.hic.p_utg.10kb.fa >> ../unitig_stats.txt

# Index & map Hi-C reads
cd ${wd}/A02-HapHiC
mkdir n2-unitigs
cd n2-unitigs
bwa index ${wd}/A01-PhasedAssembly/${species}_hifi_hi-c_n2.asm.hic.p_utg.10kb.fa
bwa mem -5SP -t 28 ${wd}/A01-PhasedAssembly/${species}_hifi_hi-c_n2.asm.hic.p_utg.10kb.fa \
    ${hic}/Zbulbosa-1A_1.fq.gz \
    ${hic}/Zbulbosa-1A_2.fq.gz \
    | samblaster | samtools view -b -@28 -S -h -F 3340 -o ${species}_unitigs_HiC.bam

# Filter BAM (MAPQ >=30 and NM <=3)
filter_bam ${species}_unitigs_HiC.bam 30 --nm 3 --threads 28 \
  | samtools view -b -@28 -o ${species}_unitigs_HiC.filtered.bam

#### Check alignment statistics with flagstat you should get ~25% of reads mapping
samtools flagstat ${species}_unitigs_HiC.filtered.bam > flagstat.report

#### Use --quick_view feature to see how well  Hi-C has mapped (the 11 is the number of chromosomes specified
#### but for quick_view this is ignored but you still need to put an arbitary value)
haphic pipeline ${wd}/A01-PhasedAssembly/${species}_hifi_hi-c_n2.asm.hic.p_utg.fa ${species}_unitigs_HiC.filtered.bam 11 --quick_view
