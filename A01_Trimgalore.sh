# =====================================================================
# Trimgalore
# Samples: paired end data with file names modified like below
# S1-R1.fastq.gz, S1-R2.fastq.gz, ..., S12-R1.fastq.gz, S12-R2.fastq.gz
# =====================================================================

#!/bin/bash
#BSUB -J "job[1-12]" trimgalore_paired
#BSUB -o output/J12_trimgalorePaired_%J-%I.out
#BSUB -e output/J12_trimgalorePaired_%J-%I.err
#BSUB -q short
#BSUB -n 1 
#BSUB -R rusage[mem=4000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu 

module load trimgalore/0.4.3
module load cutadapt/1.8.1
module load fastqc/0.11.2

trim_galore --paired --quality 20 --fastqc --stringency 1 ../S$LSB_JOBINDEX-R1.fastq.gz ../S$LSB_JOBINDEX-R2.fastq.gz --output_dir ../01_trimmed


