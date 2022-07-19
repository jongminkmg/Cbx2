#!/bin/bash
#BSUB -J CellRanger_count_S01
#BSUB -o output/countS01-%J.out
#BSUB -e output/countS01-%J.err
#BSUB -q big-multi
#BSUB -n 8
#BSUB -R rusage[mem=70000]
#BSUB -u jongminkim@molbio.mgh.harvard.edu

module load cellranger/6.0

cellranger count --id=S01 --fastqs=../Fastq --sample=S01 --transcriptome=../../../00_genome/CellRangerMouse/refdata-cellranger-mm10-3.0.0 --expect-cells=5000 --localcores=8 --localmem=64
