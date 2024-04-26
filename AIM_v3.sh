#!/bin/bash

# FILENAME: AIM_v3.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --job-name AIM_v3_array

#SBATCH --output=/home/dryals/ryals/ahb/outputs/dump.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/dump.out

#Dylan Ryals 26 OCT 2023
#last edit   13 MAR 2024

#this version reads a filename from a file to be more modular
    #this could be some sort of command line argument or named pipe, but eh

#launch with sbatch --array=1-16 AIM_v3.sh

date

module load bioinfo bcftools r

n=$( echo $SLURM_ARRAY_TASK_ID )
log=/depot/bharpur/data/projects/ryals/ahb/outputs/aim.out

cd /home/dryals/ryals/ahb/aim
filename=$( cat ref_filename.txt )
echo "reading $filename ..."

cd /scratch/bell/dryals/ahb/aim
mkdir -p chr${n}
cd chr${n}

echo "starting chr $n" >> $log

bcftools view ../../$filename -r $n -Ob -o chr${n}refs.bcf.gz
bcftools index -c chr${n}refs.bcf.gz

for pop in A C M O
#for pop in A C L M O U Y
do
    echo "starting $n $pop ..."
    bcftools view chr${n}refs.bcf.gz -S ~/ryals/ahb/references/${pop}.txt -Ou | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%AF\n' -o ${pop}.frq
    
    awk '{print $3}' ${pop}.frq > ${pop}.tmp
done

paste A.frq C.tmp M.tmp O.tmp > chr${n}.popfrq
#paste A.frq C.tmp L.tmp M.tmp O.tmp U.tmp Y.tmp > chr${n}.popfrq
rm *.tmp *.frq

echo "calculating Ia for chr $n" >> $log

Rscript --vanilla --silent /depot/bharpur/data/projects/ryals/ahb/aimIa_v2.R $n
    #v2 runtime: ~3:40
    #v1 runtime: ~9:30

#report to logs
echo "FINISHED CHR $n" >> $log
echo "    FINISHED CHR $n" >> /depot/bharpur/data/projects/ryals/ahb/outputs/pipeline.out

echo "done"
date
