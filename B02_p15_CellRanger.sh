#!/bin/bash
#BSUB -J CRanger_count_do18195
#BSUB -o output/do18195-%J.out
#BSUB -e output/do18195-%J.err
#BSUB -q big
#BSUB -n 5
#BSUB -R rusage[mem=30000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu

module add cellranger/3.0.2

cellranger count --id=P15_do18195 --fastqs=../11_fastq --sample=do18195 --transcriptome=../../../00_genome/CellRangerMouse/refdata-cellranger-mm10-3.0.0 --expect-cells=4000 --localcores=5 --localmem=29
