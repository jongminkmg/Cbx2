# ======================================================
# Sam to Bam (q10, f2)  / Picard dup / Sort / Index / QC
# ======================================================

#!/bin/bash
#BSUB -J "job[1-12]" StoB_Sort_Index
#BSUB -o output/StoB_SortIndex_-%J-%I.out
#BSUB -e output/StoB_SortIndex_-%J-%I.err
#BSUB -q medium
#BSUB -n 4
#BSUB -R rusage[mem=8000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu

module load samtools/1.4.1
module load picard-tools-2.17.1
module load java/1.8.0_181


samtools flagstat ../01_trimmed/ST$LSB_JOBINDEX.sam > ../01_trimmed/Flagstat_ST$LSB_JOBINDEX.txt

samtools view -Sb -q 10 ../01_trimmed/ST$LSB_JOBINDEX.sam > ../01_trimmed/ST$LSB_JOBINDEX-q10.bam

samtools flagstat ../01_trimmed/ST$LSB_JOBINDEX-q10.bam > ../01_trimmed/Flagstat_ST$LSB_JOBINDEX-q10.txt

samtools view -b -f 2 ../01_trimmed/ST$LSB_JOBINDEX-q10.bam > ../01_trimmed/ST$LSB_JOBINDEX-q10f2.bam

samtools flagstat ../01_trimmed/ST$LSB_JOBINDEX-q10f2.bam > ../01_trimmed/Flagstat_ST$LSB_JOBINDEX-q10f2.txt

samtools sort -o ../01_trimmed/ST$LSB_JOBINDEX-q10f2Srt.bam ../01_trimmed/ST$LSB_JOBINDEX-q10f2.bam

samtools index ../01_trimmed/ST$LSB_JOBINDEX-q10f2Srt.bam

java -Xmx4G -jar /apps/software/picard/2.6.0-Java-1.8.0_161/picard.jar MarkDuplicates REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=300000 \
I=../01_trimmed/ST$LSB_JOBINDEX-q10f2Srt.bam \
O=../01_trimmed/ST$LSB_JOBINDEX-q10f2SrtPcd.bam \
M=../01_trimmed/PcdMetrics-ST$LSB_JOBINDEX-q10f2Srt.txt

samtools index ../01_trimmed/ST$LSB_JOBINDEX-q10f2SrtPcd.bam

samtools flagstat ../01_trimmed/ST$LSB_JOBINDEX-q10f2SrtPcd.bam > ../01_trimmed/Flagstat_ST$LSB_JOBINDEX-q10f2SrtPcd.txt

# -q 10: discard reads with MAPQ smaller than 10
# -f 2: only get properly paired reads
                                    
# rm: need to delete all intermediate files after checking them. 
