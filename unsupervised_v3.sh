#!/bin/bash

# FILENAME: unsupervised_v3.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=2:00:00
#SBATCH --job-name usadmix-array
#SBATCH --partition cpu
#SBATCH --output=/home/dryals/ryals/ahb/outputs/er-usadmix.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/er-usadmix.out

#Dylan Ryals 26 OCT 2023
#last edit

#launch with: sbatch --array=4-9 unsupervised_v3.sh

date

n=$( echo $SLURM_ARRAY_TASK_ID )
log=/home/dryals/ryals/ahb/outputs/usadmix.out
#this is the latest file outputed by plink. allows this script to be a bit more modular
filename=$( cat /scratch/bell/dryals/ahb/plink/plink_admix_filename.txt )
cd /scratch/bell/dryals/ahb/admix/unsupervised

echo "running admixture $n ..." >> $log
    ADMIX=/home/dryals/bharpur/apps/admixture/admixture
    
    
    #try CV 20        
    $ADMIX --cv=20 ../../plink/${filename}.bed $n -j4 > ${filename}.${n}.out

echo "$n DONE" >> $log




#run this later to pull CV values from logs 
#grep "CV" admix*out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > cv.txt

#grep "CV" reference*out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > cvREF.txt

