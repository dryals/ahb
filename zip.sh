#!/bin/bash

# FILENAME: zip.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=1-00:00:00
#SBATCH --job-name zip
#SBATCH --output=/home/dryals/ryals/ahb/outputs/zip.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/zip.out

cd $CLUSTER_SCRATCH/ahb/sra

tar -czvf ../ryals-fikere-impute-fastq.tar *

echo "DONE"


