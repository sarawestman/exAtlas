#!/bin/bash

##SBATCH -w franklin
##SBATCH -p core
##SBATCH -t 48:00:00 
#SBATCH -p nolimit
#SBATCH -t unlimited 
#SBATCH --mem 200GB
#SBATCH -n 21
#SBATCH --mail-type FAIL
#SBATCH --mail-user=sara.westman@umu.se
#SBATCH --account u2015011

module load R
      
Rscript $@ 20 # set nCores for parallelisation
