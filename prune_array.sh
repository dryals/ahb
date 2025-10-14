#!/bin/bash

# FILENAME: prune_array.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --partition cpu
#SBATCH --job-name prune_array.sh

#SBATCH --output=/home/dryals/ryals/ahb/outputs/dump_prune.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/dump_prune.out

#Dylan Ryals 21 JAN 2025
#last edited 02 OCT 2025

#this version for ahb post revision

date
echo "---------------"
module load r

#cd $CLUSTER_SCRATCH/ahb/aim
n=$( echo $SLURM_ARRAY_TASK_ID )
log=~/ryals/ahb/outputs/prune.out

cd ~/ryals/ahb
echo "starting chr $n..." >> $log

Rscript --vanilla --silent LDprune_p.R $n
( flock -x 9 
    echo "FINISHED CHR $n" >> $log
    echo "    FINISHED CHR $n" >> ~/ryals/ahb/outputs/pipeline.out
) 9> ~/ryals/ahb/.PRUNEwritelock

echo "---------------"
date
