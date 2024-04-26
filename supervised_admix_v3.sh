#!/bin/bash

# FILENAME: supervised_admix_v3.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --time=3:00:00
#SBATCH --job-name supadmix_v3

#SBATCH --output=/home/dryals/ryals/ahb/outputs/supadmix.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/supadmix.out

#Dylan Ryals 26 OCT 2023
#

#admixture

date

cd /scratch/bell/dryals/ahb/admix/supervised
#.pop file is needed to determine populations, created with makeAdmixPop.R script
filename=$( cat /scratch/bell/dryals/ahb/plink/plink_admix_filename.txt )

    echo "running admixture..."
    
        ADMIX=/depot/bharpur/apps/admixture/admixture
        #remember to change this for the correct k value!!
        #try cv 20
        $ADMIX ../../plink/${filename}.bed 4 -j8 --cv=20 --supervised > ${filename}.out
        
date

echo "DONE"
