#!/bin/bash -l
set -eu

# vars
proj=u2015011
email=sara.westman@umu.se
in=/mnt/picea/storage/data/aspseq/nstreet/exAtlas_chemotypes/UKPROJ1/GB/GB_TR/PJ_GB/seq/auto/1100/project/X204/X204SC2304/X204SC23041202-Z01-F005/X204SC23041202-Z01-F005.RRA.PE150.20230510135650/data_release/X204SC23041202-Z01-F005/01.RawData
out=/mnt/picea/home/swestman/Sara/exAtlas/res/RNAseq
reference=/mnt/picea/storage/reference
singularity=/mnt/picea/projects/singularity/kogia
start=2
end=7
sortmerna_fasta=$reference/rRNA/sortmerna/v4.3.4/fasta/smr_v4.3_fast_db.fasta
sortmerna_inx=$reference/rRNA/sortmerna/v4.3.4/indices/smr_v4.3_fast_db
trimmomatic_adapter=$reference/Illumina/adapters/TruSeq3-PE.fa
salmon_index=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/Potra02_transcripts_with-decoy_salmon-version-1-dot-8-dot0

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# prep sample name file (to only inc. my files)
for fwd in "${in}/X*/*_1.fq.gz"; do
  echo ${fwd} > $out/RNA_samples.txt
done

# run
for f in `cat $out/RNA_samples.txt`; do
  bash ./runRNASeqPreprocessing.sh -s $start -e $end \
  -x $sortmerna_inx -X $sortmerna_fasta \
  -A $trimmomatic_adapter \
  -S $salmon_index -C \
  $proj $email $singularity $f ${f/_1.fq.gz/_2.fq.gz} $out
done

