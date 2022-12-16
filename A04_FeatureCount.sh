# ===============================================
# Featurecount to get transcript counts
# ===============================================

#!/bin/bash
#BSUB -J "job[1-12]" Feature-count 
#BSUB -o output/FCount_1-12_-%J-%I.out
#BSUB -e output/FCount_1-12_-%J-%I.err
#BSUB -q medium 
#BSUB -n 4 
#BSUB -M 8000 
#BSUB -R rusage[mem=8000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu 


module load subread/1.5.0 

# SMARTer Stranded Total RNA-Seq Kit v2 -Pico (Takara), Cbx2 23KRA - use strand '2'
# SMART-Seq mRNA LP (Takara), Cbx2 iKO - use strand '0'

# Per gene basis counts
featureCounts -T 4 -p -s 0 -t exon -g gene_id -a ../../../00_genome/mm10_refFlat_noChr.gtf -o ../02_sam/ST$LSB_JOBINDEX-FcountS0.txt ../02_sam/ST$LSB_JOBINDEX-q10f2Srt.bam

featureCounts -T 4 -p -s 2 -t exon -g gene_id -a ../../../00_genome/mm10_refFlat_noChr.gtf -o ../02_sam/ST$LSB_JOBINDEX-FcountS2.txt ../02_sam/ST$LSB_JOBINDEX-q10f2Srt.bam

# Per transcript basis counts
featureCounts -T 4 -p -s 0 -t exon -g transcript_id -a ../../../00_genome/mm10_refFlat_noChr.gtf -o ../02_sam/ST$LSB_JOBINDEX-FTxcountS0.txt ../02_sam/ST$LSB_JOBINDEX-q10f2Srt.bam

featureCounts -T 4 -p -s 2 -t exon -g transcript_id -a ../../../00_genome/mm10_refFlat_noChr.gtf -o ../02_sam/ST$LSB_JOBINDEX-FTxcountS2.txt ../02_sam/ST$LSB_JOBINDEX-q10f2Srt.bam

# -T: number of threads; -p: paired-end; -s: stranded

