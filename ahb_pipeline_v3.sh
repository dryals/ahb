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
#last edit   19 NOV 2025

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

    version=oct25

    echo "---------------------"
    echo "VERSION: $version"
    echo "---------------------"
#### VERSION ###
    
#     
# echo "---------------------"
# echo "filtering"
# echo "---------------------"
# 
#     
# echo "move samples to ahb dir..."
#     #raw imputed vcf from gencove, no QC
#     cd ~/ryals/ahb
#     #extract AHB samples 
#     awk -F',' '{print $2}' ahb_metadata.csv | tail -n +2 > sample_lists/admix3.txt
#     
#     cd $CLUSTER_SCRATCH/pipeline
# 
#     bcftools view allsamp.bcf.gz -S ~/ryals/ahb/sample_lists/admix2.txt --threads $SLURM_NTASKS -Ob -o ../ahb/ahbsamples.allsites.bcf.gz
#     cd ../ahb
#     echo "    indexing..."
#     bcftools index -c ahbsamples.allsites.bcf.gz
# 
# echo "filtering samples..."
#     cd $CLUSTER_SCRATCH/ahb
#     #run QC: declare < 0.99 probability missing, MAF>0.01, remove sites with > 10% missing genotypes
#     bcftools filter ahbsamples.allsites.bcf.gz -S . \
#     -i 'GP[:0] > 0.99 | GP[:1] > 0.99 | GP[:2] > 0.99' -Ou | \
#     bcftools view -q 0.01:minor -e 'F_MISSING>0.1' --threads $SLURM_NTASKS \
#     -Ob -o samples.filter.${version}.bcf.gz
#     
#     echo "    indexing..."
#     bcftools index -c samples.filter.${version}.bcf.gz
# 
# echo "creating sample site list..."
#     bcftools query samples.filter.${version}.bcf.gz -f'%CHROM\t%POS\n' -o plink/samples.${version}.sites
# 
# echo "determining appropriate references..."
#     echo "    pulling full reference file..."
#     cd $CLUSTER_SCRATCH/ahb
#     #no multiallelic sites, only snps, keep subset of references, no contigs
#     #rename chromosomes to "1,2,3...16"
#       #quite big, ~17M sites
#     bcftools view $refs -r $chrsLong -m2 -M2 \
#         -v snps -Ou | bcftools annotate --rename-chrs $rename \
#         --threads $SLURM_NTASKS -Ob -o allRef.bcf.gz
# 
#     bcftools index -c allRef.bcf.gz
#     
#     echo "    loading into plink..."
#     plink --bcf allRef.bcf.gz --make-bed \
#         --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# \
#         --bp-space 2000 \
#         --threads $SLURM_NTASKS \
#         --out plink/allRefThin
#     
#     echo "    running admixture..."
#     ADMIX=/home/dryals/bharpur/apps/admixture/admixture   
#     cd $CLUSTER_SCRATCH/ahb/admix/unsupervised
#     $ADMIX --cv=10 $CLUSTER_SCRATCH/ahb/plink/allRefThin.bed 7 -j${SLURM_NTASKS} > allRefThin.7.out
#     
#      #WARNING: run ahb_analysis.R to output lists of pure samples...
# 
# echo "filtering references..."
# #filter references to informative sites
#     cd $CLUSTER_SCRATCH/ahb
#     bcftools view allRef.bcf.gz \
#         -T plink/samples.${version}.sites -S ~/ryals/ahb/references/pureRefs.txt -M2 -m2 \
#         --threads $SLURM_NTASKS -Ob -o reference.filter.${version}.bcf.gz
#         
#     echo "  indexing..."
#     bcftools index -c reference.filter.${version}.bcf.gz
#     
#     #kill script if the above fails
#     if [ ! -f "reference.filter.${version}.bcf.gz.csi" ]; then
#         echo "Filters Failed!"
#         exit 1
#     fi
#    
echo "---------------------"
echo "selecting sites"
echo "---------------------"
#    
echo "launching Ia script"
    #count number of samples in each population
    cd /home/dryals/ryals/ahb/references
    wc -l ?.txt | awk '{print $1}' > refN.txt
    #reset logifle
    cd /home/dryals/ryals/ahb
    mkdir -p aim
    echo -n "" > outputs/aim.out
    #specify reference file
    echo "reference.filter.${version}.bcf.gz" > aim/ref_filename.txt
    #launch the admixture array
    sbatch --array=1-16 AIM_v3.sh
# 
# echo "merging samples and references..."
#     cd /scratch/bell/dryals/ahb
#     #merge and remove new multialleles
#     bcftools merge reference.filter.${version}.bcf.gz samples.filter.${version}.bcf.gz -m snps -Ou | \
#         bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS \
#         -Ob -o admix.${version}.bcf.gz
#     echo "  indexing..."
#     bcftools index -c admix.${version}.bcf.gz
# 
# echo "counting total sites..."
#     bcftools index -c admix.${version}.bcf.gz
#     bcftools view admix.${version}.bcf.gz | grep -vc "#" > admix.${version}.sitecount

# 
#wait for Ia calculation to finish if it hasn't already
echo "waiting for Ia results (see aim.out)..."
    cd /home/dryals/ryals/ahb
    while [ $(grep "FINISHED" outputs/aim.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
    do
        sleep 10 #wait between each check
    done
    echo ""
    
echo "compiling Ia results..."    
    cd /scratch/bell/dryals/ahb/aim
    #this will hold all the aims
    cat chr*/chr*.ia | grep -v "chr" | sort -k3 -gr > aim.ia.txt

    #TODO: verify Ia is calculated correctly, investigate NA's (equally dispersed across genome?)
    
#     #output top sites in plink format -- chr:pos
#         awk 'OFS=":" {print$1, $2}' aim.ia.txt | head -n 2000 > plink_aim.${version}.txt
    #IA greater than zero
    grep -v "NA" aim.ia.txt | awk '$3>0' | awk 'OFS=":" {print$1, $2}' \
        > plink_aim.${version}.txt
    
    count=$( wc -l aim.ia.txt | awk '{print $1}')
    echo "    Calculated Ia for $count sites"
    
echo "plink: filtering whole file for AIMs ..."    
    cd $CLUSTER_SCRATCH/ahb
    #bp filter
    plink --bcf admix.${version}.bcf.gz --make-bed \
        --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# \
        --extract aim/plink_aim.${version}.txt --threads $SLURM_NTASKS --silent \
        --out plink/topaim.${version}
    
    #LD pruning 
            #extract ld data, removing references
        echo "    calculating ld..."
        cd $CLUSTER_SCRATCH/ahb/plink
        plink --bfile topaim.${version} -r2 --ld-window 1000 --ld-window-kb 50 --ld-window-r2 0.2 \
            --remove /home/dryals/ryals/ahb/references/plink_refs.txt \
            --silent --threads $SLURM_NTASKS --out preprune
            
        #run R script to generate best set of sites by Ia
            #reset logfile
            cd ~/ryals/ahb
            echo -n "" > outputs/prune.out
            #start
            sbatch --array=1-16 prune_array.sh
            #wait
            echo "    waiting for pruning (see prune.out)..."
            while [ $(grep "FINISHED" outputs/prune.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
            do
                sleep 10 #wait between each check
            done
            
            #create full output
            echo "    compiling results..."
            cd $CLUSTER_SCRATCH/ahb/aim
            #this will hold all the aims
            cat chr*/LDremove.txt > allLDremove.txt
            count=$( wc -l allLDremove.txt | awk '{print $1}')
            echo "    marked $count sites"
            
        echo "    removing pruned sites..."
        cd $CLUSTER_SCRATCH/ahb/plink
        #create new admix file
        plink --bfile topaim.${version} --make-bed --exclude ../aim/allLDremove.txt --silent \
            --threads $SLURM_NTASKS --out admix.${version}
    
    #use this plink file basename for admix scripts
    echo "admix.${version}" > plink_admix_filename.txt
    
    
# echo "plink: pulling references..."  
#     #for unsupervised reference admix
#     cd $CLUSTER_SCRATCH/ahb/plink
#     plink --bfile admix.${version} --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
#         --keep ~/ryals/ahb/references/plink_refs.txt --threads $SLURM_NTASKS --silent --out reference.${version}
#     
    #kill script if the above fails
    if [ ! -f "admix.${version}.bed" ]; then
        echo "Plink Failed!"
        exit 1
    fi

echo "---------------------"
echo "Analysis"
echo "---------------------"
    
    
echo "starting admix..."
    cd /scratch/bell/dryals/ahb
    mkdir -p admix
    cd admix
    mkdir -p unsupervised
    mkdir -p supervised
  
    #supervised
        #create pop file
        cd /home/dryals/ryals/ahb
        R --vanilla --no-save --no-echo --silent < makeAdmixPop.R
        sleep 5
        sbatch supervised_admix_v3.sh
        
        #previous CV with unpure refs: 0.32814
        #new CV witwh pure refs: 0.33835 LAME
    
# echo "plink: generating PCA..."
#     cd $CLUSTER_SCRATCH/ahb
#     plink --bcf samples.filter.${version}.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --threads $SLURM_NTASKS --silent --maf 0.05 --pca 500 --out plink/samps.${version}
#     
#     plink --bcf admix.${version}.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --silent --threads $SLURM_NTASKS --maf 0.05 --pca 500 --out plink/all.${version}
# 
# echo "starting reference admix..."
# 
#     #WARNING: this will break if an admix file already exisit for the version
#     cd $CLUSTER_SCRATCH/ahb/admix/supervised
#     while [ ! -f "admix.${version}.4.Q" ] 
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
#     sbatch --array=4 unsupervised_v3.sh
    
echo "---------------------"
echo "---------------------"
echo "my work here is done"
date
