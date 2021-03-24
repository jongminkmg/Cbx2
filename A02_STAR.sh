# =================================================
# Mapping trimmed RNA-Seq fastq files by using STAR
# =================================================

# Note, with unknown error, STAR does not run as array job in Erisone. Hence single jobs.

#!/bin/bash
#BSUB -J S1_map_STAR
#BSUB -o output/STAR-%J-%I.out
#BSUB -e output/STAR-%J-%I.err
#BSUB -q big 
#BSUB -n 5 
#BSUB -R rusage[mem=35000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu 

module load star/2.5.3

star --runThreadN 5 --genomeDir ../../../00_genome/STAR_mm10_genome --readFilesIn ../01_trimmed/S1-R1_val_1.fq ../01_trimmed/S1-R2_val_2.fq --outFileNamePrefix ../01_trimmed/S1-


