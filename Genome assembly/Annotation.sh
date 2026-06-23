module load Anaconda3/2022.05

wd=/mnt/parscratch/users/bo1lpg/NEOF-project/sequencing/2025-assembly/
species='Zuloagaea_bulbosa'
gffread=/mnt/parscratch/users/bo1lpg/NEOF-project/sequencing/annotation/gffread

# export HELIXER_SIF=${wd}/Helixer/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif
#### Done once to fetch plant training model
# apptainer exec ${HELIXER_SIF} fetch_helixer_models.py --lineage land_plant

#### Create directory and copy assembly
mkdir ${wd}/A06-Helixer_annotation
cd ${wd}/A06-Helixer_annotation
cp ${wd}/A04-YAHS/n2-unitigs-no-correction-clean/yahs.out_scaffolds_final.fa ./${species}_assembly.fa

#### Run Helixer
apptainer exec ${HELIXER_SIF} Helixer.py --fasta-path ${species}_assembly.fa \
  --lineage land_plant --batch-size 512 --gff-output-path Zuloagaea_bulbosa.gff3

#### Get cds and aa files
${gffread}/gffread -g ${species}_assembly.fa -x ${species}.cds.fa Zuloagaea_bulbosa.gff3
${gffread}/gffread -g ${species}_assembly.fa -y ${species}.aa.fa Zuloagaea_bulbosa.gff3

#### Run BUSCO
source activate BUSCO
busco --in  ${species}.cds.fa \
      --out annotation \
      --lineage_dataset poales_odb10 \
      --mode transcriptome \
      --cpu 18 \
