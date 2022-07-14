#!/bin/bash
#BSUB -J CellRanger_count
#BSUB -o output/test-%J.out
#BSUB -e output/test-%J.err
#BSUB -q big-multi
#BSUB -n 8
#BSUB -R rusage[mem=30000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu

module add cellranger/2.0.0

cellranger count --id=S2WTmlti --fastqs=../97_fastq_S2 --sample=JOc4 --transcriptome=../../../00_genome/CellRangerMouse/refdata-cellranger-mm10-1.2.0 --expect-cells=3000 --localcores=8 --localmem=29
