#!/bin/bash

# FILENAME: 1_array_filter_v3.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=3:00:00
#SBATCH --job-name array_filter_v3
#SBATCH --partition cpu
#SBATCH --output=/home/dryals/ryals/ahb/outputs/null.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/null.out

#Dylan Ryals 26 OCT 23

#launch with: sbatch --array=1-2 1_array_filter_v3.sh

#this will filter input files using a job array
#should take <1hr

module load bioinfo bcftools

#this variable is the ID for the job array
n=$( echo $SLURM_ARRAY_TASK_ID )
samp=$( sed -n ${n}p sample_lists/vcf_files_needed.txt )
log=/home/dryals/ryals/ahb/outputs/array.out

echo "working on $samp" >> $log
echo "  working on $samp" >> /home/dryals/ryals/ahb/outputs/grand.out
echo "filtering sample $n ..." >> $log

    chrs=$( cat /home/dryals/ryals/diversity/chrfilter.txt | tr '\n' ',' )
    rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
    
    cd /scratch/bell/dryals/ahb
    
    #keep no contigs, only bialelic, rename chrs
    bcftools view $samp -M2 -v snps -f .,PASS -r $chrs -Ou | bcftools norm -m +snps -Ou | bcftools view $samp -M2 -v snps -Ou | bcftools annotate --rename-chrs $rename --threads 4 -Ob -o sample$n.bcf.gz
    
    bcftools index -c sample$n.bcf.gz
    
    #consider removing sites with low GP
    
#ouptut to various logs 
echo "DONE with $n" >> $log
echo "    DONE with $n" >> /home/dryals/ryals/ahb/outputs/grand.out

#report to sentinel file so GRAND can stop waiting
echo "DONE with $n" >> /home/dryals/ryals/ahb/sentinel.txt


