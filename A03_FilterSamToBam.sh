# ==============================================
# Sam to Bam (filter w/ q10, f2)  / Sort / Index
# ==============================================

#!/bin/bash
#BSUB -J "job[1-12]" StoB_Sort_Index
#BSUB -o output/StoB_SortIndex_-%J-%I.out
#BSUB -e output/StoB_SortIndex_-%J-%I.err
#BSUB -q medium 
#BSUB -n 1 
#BSUB -R rusage[mem=8000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu 

# Note skip Pcard duplicate removal since this is RNA-seq
# Run after moving sam files into ../02_sam folder from ../01_trimmed

module load samtools/1.4.1


samtools flagstat ../02_sam/S$LSB_JOBINDEX-Aligned.out.sam > ../02_sam/Flagstat_S$LSB_JOBINDEX.txt

# Mapquality filter - obtain uniquely mapped reads
# I believe STAR map quality score 255 = 1 map, 3 = 2 map, 2 = 3 map, ...
# MAPQ = int(-10*log10(1-1/Nmap)) 
samtools view -Sb -q 10 ../02_sam/S$LSB_JOBINDEX-Aligned.out.sam > ../02_sam/ST$LSB_JOBINDEX-q10.bam

samtools flagstat ../02_sam/ST$LSB_JOBINDEX-q10.bam > ../02_sam/Flagstat_ST$LSB_JOBINDEX-q10.txt

# Filter to obtain proper pairs
samtools view -b -f 2 ../02_sam/ST$LSB_JOBINDEX-q10.bam > ../02_sam/ST$LSB_JOBINDEX-q10f2.bam

samtools flagstat ../02_sam/ST$LSB_JOBINDEX-q10f2.bam > ../02_sam/Flagstat_ST$LSB_JOBINDEX-q10f2.txt

samtools sort -o ../02_sam/ST$LSB_JOBINDEX-q10f2Srt.bam ../02_sam/ST$LSB_JOBINDEX-q10f2.bam

samtools index ../02_sam/ST$LSB_JOBINDEX-q10f2Srt.bam

