#!/bin/bash
set -eu

in=/mnt/picea/home/swestman/Sara/exAtlas/res/RNAseq/fastqc/trimmomatic
out=/mnt/picea/home/swestman/Sara/exAtlas/res/RNAseq/multiqc/trimmomatic
singularity=/mnt/picea/projects/singularity/kogia/multiqc_1.11.sif

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# run
sbatch ./runMultiQC.sh $singularity $in $out