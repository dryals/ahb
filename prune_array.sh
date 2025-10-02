#!/bin/bash

# FILENAME: prune_array.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --partition cpu
#SBATCH --job-name prune_array.sh

#SBATCH --output=/home/dryals/ryals/admixPipeline/outputs/dump_prune.out
#SBATCH --error=/home/dryals/ryals/admixPipeline/outputs/dump_prune.out

#Dylan Ryals 21 JAN 2025
#last edited 02 OCT 2025

#this version for ahb post revision

date
echo "---------------"
module load r

cd $CLUSTER_SCRATCH/pipeline/aim
n=$( echo $SLURM_ARRAY_TASK_ID )
log=/depot/bharpur/data/projects/ryals/admixPipeline/outputs/prune.out

cd ~/ryals/admixPipeline
echo "starting chr $n..." >> $log

Rscript --vanilla --silent LDprune_p.R $n

echo "FINISHED CHR $n" >> $log

echo "---------------"
date
