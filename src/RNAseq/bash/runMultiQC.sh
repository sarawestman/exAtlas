#!/bin/bash
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mem=16GB
#SBATCH --mail-type=END,FAIL
#SBATCH -t 12:00:00
#SBATCH --mail-user=sara.westman@umu.se
#SBATCH --account u2015011

## stop on error but be verbose
set -eux

# source helpers
source ./functions.sh

USAGETXT=\
"
Usage: $(basename $0) <singularity container> <analysis directory> <output directory>
"

## arguments
[[ $# -ne 3 ]] && abort "This script takes three arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be an existing mutiqc singularity container file."

## enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

[[ ! -d $2 ]] && abort "The second argument needs to be an existing analysis directory."

[[ ! -d $3 ]] && abort "The third argument needs to be an existing output directory."

## start
singularity exec $1 multiqc -o $3 $2
