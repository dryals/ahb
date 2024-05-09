#!/bin/bash

# FILENAME: sra_upload.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name sra
#SBATCH --output=/home/dryals/ryals/ahb/outputs/sra_upload.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/sra_upload.out

echo "moving..."
cd /scratch/bell/dryals/ahb/sra_upload

echo "lftp..."

lftp -u subftp,SniappegEtnurak3 ftp-private.ncbi.nlm.nih.gov << EOF

cd uploads/dylan.k.ryals_gmail.com_q46tAQoz/impute
mput *

bye
EOF


echo "DONE"


