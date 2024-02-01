#!/bin/bash

#SBATCH -w franklin
#SBATCH -p core
#SBATCH -t 48:00:00 
#SBATCH --mem 10GB
#SBATCH -n 5
#SBATCH --mail-type FAIL
#SBATCH --mail-user=sara.westman@umu.se
#SBATCH --account u2015011

~/Git/exAtlas/src/sirius/bin/sirius -i ~/Git/exAtlas/src/doc/Neg_QC11_MSMS_10_20_40_Iter1.mgf -o ~/Git/exAtlas/src/doc/my_test5 formula --elements-enforced C,H,O -p qtof fingerprint structure --database ALL_BUT_INSILICO compound-classes write-summaries --output ~/Git/exAtlas/src/doc/my_sum5

