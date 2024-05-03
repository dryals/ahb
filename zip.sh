#!/bin/bash

# FILENAME: zip.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --job-name zip
#SBATCH --output=/home/dryals/ryals/ahb/outputs/zip.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/zip.out

cd $CLUSTER_SCRATCH/ahb

tar -czvf ryals-fikere-impute-fastq.tar sra

echo "DONE"


