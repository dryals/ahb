#!/bin/bash

# FILENAME: ahb_pipeline_v2.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --time=1-00:00:00
#SBATCH --partition cpu
#SBATCH --job-name ahb_pipeline_v2.sh
#SBATCH --output=/home/dryals/ryals/ahb/outputs/pipeline.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/pipeline.out

#Dylan Ryals 15 MAR 2024
#last edit   30 SEP 2025

#added file versioning 
#revisions in 2025 

date

module load biocontainers bcftools plink r

#set some global vars for later use 
    chrsLong=$( cat /home/dryals/ryals/diversity/chrfilter.txt | tr '\n' ',' )
    refs=/depot/bharpur/data/popgenomes/HarpurPNAS/output_snp.vcf.gz
    rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )
    
#### VERSION ###
    #version will be appended to all filenames to keep versions straight

    version=sep25

    echo "---------------------"
    echo "VERSION: $version"
    echo "---------------------"
#### VERSION ###
    
    
echo "---------------------"
echo "filtering"
echo "---------------------"

    
echo "move samples to ahb dir..."
    #raw imputed vcf from gencove, no QC
    cd ~/ryals/ahb
    #extract AHB samples 
    awk -F',' '{print $2}' ahb_metadata.csv | tail -n +2 > sample_lists/admix3.txt
    
    cd $CLUSTER_SCRATCH/pipeline

    bcftools view allsamp.bcf.gz -S ~/ryals/ahb/sample_lists/admix2.txt --threads $SLURM_NTASKS -Ob -o ../ahb/ahbsamples.allsites.bcf.gz
    cd ../ahb
    echo "    indexing..."
    bcftools index -c ahbsamples.allsites.bcf.gz
# 
# echo "filtering samples..."
#     cd $CLUSTER_SCRATCH/ahb
#     #run QC: declare < 0.99 probability missing, MAF>0.01, remove sites with > 10% missing genotypes
#     bcftools filter ahbsamples.allsites.bcf.gz -S . -i 'GP[:0] > 0.99 | GP[:1] > 0.99 | GP[:2] > 0.99' -Ou | bcftools view -q 0.1:minor -e 'F_MISSING>0.1' --threads $SLURM_NTASKS -Ob -o ../ahb/samples.filter.${version}.bcf.gz
#     
#     echo "    indexing..."
#     bcftools index -c samples.filter.${version}.bcf.gz
#     
# echo "creating site list..."
#     bcftools query samples.filter.${version}.bcf.gz -f'%CHROM\t%POS\n' -o plink/samples.${version}.sites

# echo "pulling reference..."
#     cd $CLUSTER_SCRATCH/ahb
#     #no multiallelic sites, only snps, keep subset of references, no contigs, rename chromosomes to "1,2,3...16"
#     bcftools view $refs -S /home/dryals/ryals/ahb/references/all_refs.txt -r $chrsLong -M2 -v snps -Ou | bcftools annotate --rename-chrs $rename --threads $SLURM_NTASKS -Ob -o reference.bcf.gz
# 
#     bcftools index -c reference.bcf.gz
# 
# echo "filtering references..."
# #filter references to informative sites
#     cd $CLUSTER_SCRATCH/ahb
#     bcftools view reference.bcf.gz -T plink/samples.${version}.sites --threads $SLURM_NTASKS -Ob -o reference.filter.${version}.bcf.gz
#     echo "  indexing..."
#     bcftools index -c reference.filter.${version}.bcf.gz
#     
#     #kill script if the above fails
#     if [ ! -f "reference.filter.${version}.bcf.gz.csi" ]; then
#         echo "Filters Failed!"
#         exit 1
#     fi
#    
# echo "---------------------"
# echo "selecting sites"
# echo "---------------------"
#    
#    
#    
# echo "launching Ia script"
#     #count number of samples in each population
#     cd /home/dryals/ryals/ahb/references
#     wc -l ?.txt | awk '{print $1}' > refN.txt
#     #reset logifle
#     cd /home/dryals/ryals/ahb
#     mkdir -p aim
#     echo -n "" > outputs/aim.out
#     #specify reference file
#     echo "reference.filter.${version}.bcf.gz" > aim/ref_filename.txt
#     #launch the admixture array
#     sbatch --array=1-16 AIM_v3.sh
# 
# echo "merging samples and references..."
#     cd /scratch/bell/dryals/ahb
#     #merge and remove new multialleles
#     bcftools merge reference.filter.${version}.bcf.gz samples.filter.${version}.bcf.gz -m snps -Ou | bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS -Ob -o admix.${version}.bcf.gz
#     echo "  indexing..."
#     bcftools index -c admix.${version}.bcf.gz
# 
# #wait for Ia calculation to finish if it hasn't already
# echo "waiting for Ia results (see aim.out)..."
#     #WARNING: this fails if ia runs to quickly, need to write-protect an output file...
#     cd /home/dryals/ryals/ahb
#     while [ $(grep "FINISHED" outputs/aim.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
#     do
#         sleep 10 #wait between each check
#     done
#     echo ""
#     
# echo "compiling Ia results..."    
#     cd /scratch/bell/dryals/ahb/aim
#     #this will hold all the aims
#     cat chr*/chr*.ia | grep -v "chr" | sort -k3 -gr > aim.${version}.txt
#    
#     #output top sites in plink format -- chr:pos
#     awk 'OFS=":" {print$1, $2}' aim.${version}.txt | head -n 5000 > plink_aim.${version}.txt
#     count=$( wc -l aim.${version}.txt | awk '{print $1}')
#     echo "    Calculated Ia for $count sites"
#     
# echo "plink: filtering whole file for AIMs ..."    
#     cd /scratch/bell/dryals/ahb
#     #bp filter
#     plink --bcf admix.${version}.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --extract aim/plink_aim.${version}.txt --bp-space 500 --threads $SLURM_NTASKS --silent --out plink/admix.${version}
#     
# echo "plink: pulling references..."  
#     #for unsupervised reference admix
#     cd $CLUSTER_SCRATCH/ahb/plink
#     plink --bfile admix.${version} --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --keep ~/ryals/ahb/references/plink_refs.txt --threads $SLURM_NTASKS --silent --out reference.${version}
#     
#     #use this plink file basename for admix scripts
#     echo "admix.${version}" > plink_admix_filename.txt
#     #echo "sampadmix" > plink_admix_filename.txt
#     
#     #kill script if the above fails
#     if [ ! -f "admix.${version}.bed" ]; then
#         echo "Plink Failed!"
#         exit 1
#     fi
# 
# echo "---------------------"
# echo "Analysis"
# echo "---------------------"
#     
#     
# echo "starting admix..."
#     cd /scratch/bell/dryals/ahb
#     mkdir -p admix
#     cd admix
#     mkdir -p unsupervised
#     mkdir -p supervised
#   
#     #supervised
#         #create pop file
#         cd /home/dryals/ryals/ahb
#         R --vanilla --no-save --no-echo --silent < makeAdmixPop.R
#         sleep 5
#         sbatch supervised_admix_v3.sh
#     
# echo "plink: generating PCA..."
#     cd $CLUSTER_SCRATCH/ahb
#     plink --bcf samples.filter.${version}.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --threads $SLURM_NTASKS --silent --maf 0.05 --pca 500 --out plink/samps.${version}
#     
#     plink --bcf admix.${version}.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --silent --threads $SLURM_NTASKS --maf 0.05 --pca 500 --out plink/all.${version}
#     
# 
# echo "starting reference admix..."
# 
#     #WARNING: this will break if an admix file already exisit for the version
#     #wait for admix to finish
#     cd $CLUSTER_SCRATCH/admix/supervised
#     while [ ! -f "admix.${version}.4.Q" ] \
#     do
#         sleep 10 #wait between each check
#     done
#     
#     #overwrite baseneame
#     cd $CLUSTER_SCRATCH/ahb/plink
#     echo "reference.${version}" > plink_admix_filename.txt
#     #reset log file
#     cd /home/dryals/ryals/ahb
#     echo -n "" > outputs/usadmix.out
#     #launch the admixture array
#     sbatch --array=2-9 unsupervised_v3.sh
    
#ending output
echo "---------------------"
echo "---------------------"
echo "my work here is done"
date
