#!/bin/bash

##SBATCH -w franklin
#SBATCH -p core
#SBATCH -t 48:00:00 
#SBATCH --mem 5GB
#SBATCH -n 1
#SBATCH --mail-type FAIL
#SBATCH --mail-user=sara.westman@umu.se
#SBATCH --account u2015011

# https://hub.docker.com/r/wishartlab/cfmid - info

MY_SMILE=${1}
MY_MODE=${2}
OUT_FILE=${3}

singularity exec -B /mnt:/mnt ~/Git/exAtlas/src/cfm/cfm_img.sif cfm-predict ${MY_SMILE} 0.001 /trained_models_cfmid4.0/${MY_MODE}/param_output.log /trained_models_cfmid4.0/${MY_MODE}/param_config.txt 1 ${OUT_FILE}
