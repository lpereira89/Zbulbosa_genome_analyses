module load Anaconda3/2022.05
source activate assembly

#### Directories and files
wd=/mnt/parscratch/users/bo1lpg/NEOF-project/sequencing/2025-assembly
hifi=${wd}/raw/hifi
hic=${wd}/raw/hic
assembly=${wd}/A01-PhasedAssembly/Zuloagaea_bulbosa_hifi_hi-c_n2.asm.hic.p_utg.fa 
index=${wd}/A01-PhasedAssembly/Zuloagaea_bulbosa_hifi_hi-c_n2.asm.hic.p_utg.fa.fai 
species='Zuloagaea_bulbosa'
juicer_dir=${wd}/HapHiC/utils

#### Create directories
mkdir -p ${wd}/A04-YAHS/n2-unitigs-no-correction-clean
cd ${wd}/A04-YAHS/n2-unitigs-no-correction-clean

#### Map HiC reads to the assembly
bwa index ${assembly}
bwa mem -5SP -t 28 ${assembly} \
    ${hic}/Zbulbosa-1A_1_clean.fq.gz \
    ${hic}/Zbulbosa-1A_2_clean.fq.gz \
    | samblaster | samtools view -b -@28 -S -h -F 3340 -o ${species}_unitigs_HiC.bam

# Filter BAM (MAPQ >=30 and NM <=3)
filter_bam ${species}_unitigs_HiC.bam 30 --nm 3 --threads 28 \
    | samtools view -b -@28 -o ${species}_unitigs_HiC.filtered.bam
source deactivate

#### Run YAHS pipeline
source activate yahs
yahs ${assembly} ${species}_unitigs_HiC.filtered.bam --no-contig-ec
source deactivate

#### Prepare visualisation files
module load SAMtools/1.16.1-GCC-11.3.0
module load Java/17.0.4

(${juicer_dir}/juicer pre yahs.out.bin yahs.out_scaffolds_final.agp ${index} | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) \
  && (mv alignments_sorted.txt.part alignments_sorted.txt)
(java -jar -Xmx32G ${juicer_dir}/juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic.part yahs.out_scaffolds_final.chrom.sizes) \
  && (mv out.hic.part out.hic)

${juicer_dir}/juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp ${index} >out_JBAT.log 2>&1
(java -jar -Xmx32G ${juicer_dir}/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part \
  <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

#### After manual curation, to generate the curated assembly file
# juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp contigs.fa
