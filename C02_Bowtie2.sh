
# ====================================================
# Map trimmed fastqs with Bowtie2
# ====================================================


#!/bin/bash
#BSUB -J "job[1-12]" Bowtie2
#BSUB -o output/bwt2-%J-%I.out
#BSUB -e output/bwt2-%J-%I.err
#BSUB -q normal
#BSUB -n 6
#BSUB -R rusage[mem=8000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu


module load bowtie2/2.3.1

bowtie2 --end-to-end  -p 6 -N 1 -t -x /data/kingston/jongmin/00_genome/Bowtie2_mm10_genome/mm10 -1 ../01_trimmed/S$LSB_JOBINDEX-R1_val_1.fq.gz -2 ../01_trimmed/S$LSB_JOBINDEX-R2_val_2.fq.gz -S ../01_trimmed/ST$LSB_JOBINDEX.sam

#end-to-end: 'untrimmed' or 'unclipped' alignment, so not doing 'local' alignment.
#-p: number of threads
#-N: number of mismatch in the seed
#default: map multiple locations and report only the best one / vs. -k mode: search for one mor more alignment, report each / vs. -a mode: search for and report all alignments
#-t: Print the wall-clock time
#-x: index location
