#!/bin/bash

# FILENAME: unimputed.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --job-name unimputed.sh

#SBATCH --output=/home/dryals/ryals/ahb/outputs/unimputed.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/unimputed.out

#Dylan Ryals 13 MAR 2024
#last edit   

date

module load bioinfo bcftools plink r

#set some global vars for later use 
    chrsLong=$( cat /home/dryals/ryals/diversity/chrfilter.txt | tr '\n' ',' )
    refs=/depot/bharpur/data/popgenomes/HarpurPNAS/output_snp.vcf.gz
    rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )

cd $CLUSTER_SCRATCH/ahb
    
#just pull sites with sufficient evidence to call a genotype
echo "pulling unimputed sites..."
     bcftools filter samples.filter.bcf.gz -S . -i '( FMT/DP > 3 & FMT/RC == 0 ) | ( FMT/DP > 3 & FMT/AC == 0 ) | ( FMT/AC > 0 & FMT/RC > 0)' -Ou | bcftools view -q 0.01:minor -e 'F_MISSING>0.1' --threads $SLURM_NTASKS -Ob -o samples.unimpt.bcf.gz
    echo "    indexing..." 
     bcftools index -c samples.unimpt.bcf.gz
    
    
echo "creating site list..."
    bcftools query samples.unimpt.bcf.gz -f'%CHROM\t%POS\n' -o plink/unimptSites.txt
    bcftools query samples.unimpt.bcf.gz -f'%CHROM\t%POS\t%FILTER\n' -o unimpt.sites


echo "filtering references..."
#filter references to informative sites
    bcftools view reference.bcf.gz -T plink/unimptSites.txt --threads $SLURM_NTASKS -Ob -o reference-filter-unimpt.bcf.gz
    echo "  indexing..."
    bcftools index -c reference-filter-unimpt.bcf.gz
    
    #kill script if the above fails
    if [ ! -f "reference-filter-unimpt.bcf.gz.csi" ]; then
        echo "Filters Failed!"
        exit 1
    fi
    

#WARNING: This will over-write the .ia list
    
echo "launching Ia script"
    #count number of samples in each population
    cd /home/dryals/ryals/ahb/references
    wc -l ?.txt | awk '{print $1}' > refN.txt
    #reset logifle
    cd /home/dryals/ryals/ahb
    mkdir -p aim
    echo -n "" > outputs/aim.out
    #specify reference file
    echo "reference-filter-unimpt.bcf.gz" > aim/ref_filename.txt
    #launch the admixture array
    sbatch --array=1-16 AIM_v3.sh

    
echo "merging samples and references..."
    cd /scratch/bell/dryals/ahb
    #merge and remove new multialleles
    bcftools merge reference-filter-unimpt.bcf.gz samples.unimpt.bcf.gz -m snps -Ou | bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS -Ob -o admix-unimpt.bcf.gz
    echo "  indexing..."
    bcftools index -c admix-unimpt.bcf.gz

##Ia calc no longer takes a long time, script just for reference
## echo "waiting for Ia results (see aim.out)..."
##     cd /home/dryals/ryals/ahb
##     while [ $(grep "FINISHED" outputs/aim.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
##     do
##         sleep 10 #wait between each check
##     done
#    
echo "compiling Ia results..."    
    cd /scratch/bell/dryals/ahb/aim
    #this will hold all the aims
    cat chr*/chr*.ia | grep -v "chr" | sort -k3 -gr > aim.ia.txt
    
        #rename chrs (in case I need this line again...)
        #rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
        #while read a b; do sed -i "s/$a/$b/g" aim.txt ; done < $rename
   
   
    #output top 50k sites in plink format -- chr:pos
    awk 'OFS=":" {print$1, $2}' aim.ia.txt | head -n 5000 > plink_aim-unimpt.ia.txt
    count=$( wc -l aim.ia.txt | awk '{print $1}')
    echo "Calculated Ia for $count sites"
    
echo "plink: filter whole file for AIMs ..."    
    #aims, maf, and geno
    cd /scratch/bell/dryals/ahb
    #this step should be removed when the merged file only includes sites that pass the sample MAF filter... 
    plink --bcf admix-unimpt.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --extract aim/plink_aim-unimpt.ia.txt --threads $SLURM_NTASKS --silent --out plink/admix-unimpt
    
#     #filter just samples for IA for unsupervised admix
#     plink --bcf samples.filter.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --extract aim/plink_aim.ia.txt --bp-space 1000 --threads $SLURM_NTASKS --silent --out plink/samp
    
    #use this plink file basename for admix scripts
    cd plink
    echo "admix-unimpt" > plink_admix_filename.txt
    #echo "sampadmix" > plink_admix_filename.txt
    
    #kill script if the above fails
    if [ ! -f "admix-unimpt.bed" ]; then
        echo "Plink Failed!"
        exit 1
    fi
    
echo "starting admix..."
    cd /scratch/bell/dryals/ahb
    mkdir -p admix
    cd admix
    mkdir -p unsupervised
    mkdir -p supervised
    #unsupervised
        #reset log file
        cd /home/dryals/ryals/ahb
        echo -n "" > outputs/usadmix.out
        #launch the admixture array
        sbatch --array=2-9 unsupervised_v3.sh
        
    #supervised
        #create pop file
        cd /home/dryals/ryals/ahb
        R --vanilla --no-save --no-echo --silent < makeAdmixPop.R
        sleep 5
        sbatch supervised_admix_v3.sh
    
echo "plink: generate PCA..."
    cd $CLUSTER_SCRATCH/ahb
    
    plink --bcf samples.unimpt.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --threads $SLURM_NTASKS --silent --maf 0.05 --pca 500 --out plink/samps-unimpt
    
    plink --bcf admix-unimpt.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --silent --threads $SLURM_NTASKS --maf 0.05 --pca 500 --out plink/all-unimpt
    
 
echo "my work here is done"
date
