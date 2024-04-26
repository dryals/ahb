#!/bin/bash

# FILENAME: n_mitotype.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --job-name mitotype
#SBATCH --output=/home/dryals/ryals/ahb/outputs/nmito.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/nmito.out

#Dylan Ryals 31 OCT 2023
#last edited 07 NOV 2023

#launch with:
    #sbatch --array=1-23 n_mitotype.sh
    #each array member will run 12 mitotypes
    
#this version for negishi

date

module load biocontainers bwa atram spades abyss velvet blast


#setup
    k=1 #number of mitotypes per job
    nstart=$((( $SLURM_ARRAY_TASK_ID - 1 ) * $k + 1))
    nend=$(( $nstart + $k - 1 ))
    echo "$nstart $nend"
    
    hook=~/ryals/ahb/AHB_markers.fasta
    log=~/ryals/ahb/outputs/array_mito.out
    
    #sample=$( tail -n+2 AHB_meta3.csv | sed -n "${nstart},${nend};d" | awk -F',' '{print $3}' )
    #sample=$( sed "${n}q;d" sample_lists/scratch.txt )
    #sample=$( sed "${n}q;d" sample_lists/failed_samples_2.txt )
    #sample=A3
    #k=200c1734-e9e1-4677-90ef-562c358c9d0b
    
    
#loop
echo "running mitotypes: $nstart - $nend" >> $log

#tail -n+2 AHB_meta3.csv | sed -n "${nstart},${nend} p" | awk -F',' '{print $3}' | while read sample
#sed -n "${nstart},${nend} p" sample_lists/missing.txt | while read sample
sed -n "${nstart},${nend} p" pope/popesamp.txt | while read sample
do
        cd ${CLUSTER_SCRATCH}/ahb/fastq/${sample} || continue
        k=$( ls )
        cd $k || exit
            mkdir -p lib
            mkdir -p tmp
            mkdir -p out
        #empty dirs if there was some previous run
            rm -rf tmp/* out/* lib/* *log

    #atram
    echo "    $sample building library (atram)..." >> $log
        #create libraries
        atram_preprocessor.py --cpus $SLURM_NTASKS -b ./lib/${sample} --end-1 ${k}_R1.f*q.gz --end-2 ${k}_R2.f*q.gz -t ./tmp --gzip -s 50

        #run atram
            #spades seems to work on negishi while abyss works on bell...
            #hopefully there aren't major differences!
    echo "    $sample creating contigs (atram)..." >> $log
        atram.py -b ./lib/${sample} -Q $hook -i 1 --cpus $SLURM_NTASKS -a spades --length 100 --fraction 1.0 --log-file ${sample}log -o ./out/${sample} -t ./tmp --bit-score 100
        
    #stop if atram fails...
        if [ $(ls out | wc -l | awk '{print $1}') -lt 1 ]; then
            echo "ERR : $sample aTRAM_Failed!" >> $log
            exit 1
        fi

    #blast against baits
    echo "    $sample aligning contigs to hooks (blast)..." >> $log
        makeblastdb -in $hook -dbtype nucl -parse_seqids -out lib/blast_db

        for pop in C1 A1 M3 O5
        do
            blastn -db lib/blast_db -query out/*${pop}*filter* -outfmt 6 -out out/${pop}_blast.out
            #sort -k 12 -g -r out/${pop}_blast.out | head -n 2 >> besthits.txt
        done
    #stop if blast fails...
        if [ $(ls out/*_blast.out | wc -l | awk '{print $1}') -lt 1 ]; then
            echo "ERR : $sample BLAST Failed!" >> $log
            exit 1
        fi
    #output unique contigs if all goes well.
        ( flock -x 9 
        
            echo "sampleID:${sample}" >> ~/ryals/ahb/pope/mitotype.out
            cat out/*_blast.out | sort -k 12 -g -r | awk '!a[$12]++' | head -n 10  >> ~/ryals/ahb/pope/mitotype.out
        
        ) 9> ~/ryals/ahb/.writelock
        
        echo "    $sample FINISHED" >> $log
        #cleaning
        rm -rf tmp/* #lib/*
done

echo "FINISHED mitotypes: $nstart - $nend" >> $log

date

#temporary: read all outfiles


# cat popesamp.txt | while read sample
# do
#     #sample=CA23CY04
#     echo $sample
#     id=$( grep "$sample" /home/dryals/ryals/ahb/gencove/apr_submission/sampleNames.txt | awk '{print $2}' )
#     echo -n "${id}," >> unarch.txt
# done
# 
# samples=$( cat unarch.txt )
# projid=$( cat /home/dryals/ryals/ahb/gencove/apr_submission/projID.txt ) 
# gencove projects restore-samples $projid --sample-ids $samples
# 




# cat sample_lists/special_samples.txt | while read sample
# do
#     echo "sampleID:${sample}" >> mitotype3.out
#     cat ${CLUSTER_SCRATCH}/ahb/fastq/${sample}/*/out/*_blast.out >> mitotype3.out
# done






    

