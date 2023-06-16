#!/bin/sh
#SBATCH -J RNASeq-ref-Prep
#SBATCH --partition batch
#SBATCH --ntasks=1                    # Run a single task    
#SBATCH --cpus-per-task=28             # Number of CPU cores per task
#SBATCH --time=8:00:00
#SBATCH --mem=30gb
#SBATCH --mail-user=your_email@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
date
cd $SLURM_SUBMIT_DIR
cp ./genome/g717v5_.gff ./gene.gff # copy annotation file
cp ./genome/g717v5_.fa ./genome.fa # copy genome file
module load gffread
gffread gene.gff -T -o gene.gtf # convert gff to gtf
ml STAR # load STAR and build index
STAR \ 
--runThreadN 8 \
--runMode genomeGenerate \
--genomeSAindexNbases 13 \
--genomeDir . \
--genomeFastaFiles genome.fa \
--sjdbGTFfile gene.gtf

date
