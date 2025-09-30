#!/bin/bash

# FILENAME: n_mitotype.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --partition cpu
#SBATCH --job-name mitotype
#SBATCH --output=/home/dryals/ryals/ahb/outputs/nmito.out
#SBATCH --error=/home/dryals/ryals/ahb/outputs/nmito.out

#Dylan Ryals 31 OCT 2023
#last edited 26 AUG 2024

#launch with:
    #sbatch --array=1-n n_mitotype.sh
    #where k is the number of mitotypes per job
    #each array member will run 12 mitotypes
    
#this is an array job. It will create as many identical tasks as you give it in the above --array flag
#each 'task' takes 4 cores so it has the neccesary RAM to run
#each task will attempt to allign k mitotypes from a list of sample ID's
    #you can play with these values, but k = 4-10 works well for me
    #try not to use all the cores on bharpur, leave some for others
#this only runs on negishi. If you need to run on bell, align wtih abyss instead of spades
    #and change the 'module load' command arguments

date

module load biocontainers bwa atram spades abyss velvet blast


#setup
    #this helps manage parallel processing
        #each job in the array works on k mitotypes
        #nstart and nend are the first and last mitotypes in the list a given array will process 
    k=1 #number of mitotypes per job
    nstart=$((( $SLURM_ARRAY_TASK_ID - 1 ) * $k + 1))
    nend=$(( $nstart + $k - 1 ))
    echo "$nstart $nend"
    
    #output directory for files you need to save
        #aka not scratch
        #change this to your output dir, not mine
    outdir=~/ryals/ahb
    
    #fasta file with all the haplotypes: A, C, M, O
        #you can copy this file to your own outdir
    hook=~/ryals/ahb/AHB_markers.fasta
    
    #logfile that all the arrays will report to
        #the output files in the header above will be really chaotic and not so useful
        #again, change those (header abover) and this (below) to your own outdir
    log=~/ryals/ahb/outputs/array_mito.out
    
    
#main loop
    echo "running mitotypes: $nstart - $nend" >> $log
    #sample file
        #just a list of sample ID's, aka gencove ID's
        #take just the mitotypes this task needs: nstart to nend
        #these are piped into the while loop
        #the loop uses the $sample variable to select the sample to work on
    tail -n+2 AHB_meta3.csv | sed -n "${nstart},${nend} p" | awk -F',' '{print $3}' | while read sample
    do
            #see if a fastq exists for this value of $sample, else continue to the next
            cd ${CLUSTER_SCRATCH}/ahb/fastq/${sample} || continue
            #here I confusingly overwite the $k varible oops
            #this is the long-form id for each fastq given by gencove
            k=$( ls )
            cd $k || exit
                #make temporary directories for atram
                mkdir -p lib
                mkdir -p tmp
                mkdir -p out
                #empty these dirs if there was some previous run
                rm -rf tmp/* out/* lib/* *log

        #ATRAM: build DB
        echo "    $sample building library (atram)..." >> $log
            #create blast libraries
            atram_preprocessor.py --cpus $SLURM_NTASKS -b ./lib/${sample} \
                --end-1 ${k}_R1.f*q.gz --end-2 ${k}_R2.f*q.gz -t ./tmp --gzip -s 50

        #ATRAM: create contigs
                #spades seems to work on negishi while abyss works on bell...
                #hopefully there aren't major differences!
        echo "    $sample creating contigs (atram)..." >> $log
            atram.py -b ./lib/${sample} -Q $hook -i 1 --cpus $SLURM_NTASKS -a spades \
                --length 100 --fraction 1.0 --log-file ${sample}log -o ./out/${sample} \
                -t ./tmp --bit-score 100
            
        #stop and exit if atram fails...
            if [ $(ls out | wc -l | awk '{print $1}') -lt 1 ]; then
                echo "ERR : $sample aTRAM_Failed!" >> $log
                exit 1
            fi

        #BLAST against the supplied fasta aka $hook
        echo "    $sample aligning contigs to hooks (blast)..." >> $log
            makeblastdb -in $hook -dbtype nucl -parse_seqids -out lib/blast_db

            #loop through the haplotypes and align to each of them
            for pop in C1 A1 M3 O5
            do
                blastn -db lib/blast_db -query out/*${pop}*filter* -outfmt 6 -out out/${pop}_blast.out
                #I can't remember what this did, but it's commented out
                #sort -k 12 -g -r out/${pop}_blast.out | head -n 2 >> besthits.txt
            done
            
            #stop and exit if blast fails...
            if [ $(ls out/*_blast.out | wc -l | awk '{print $1}') -lt 1 ]; then
                echo "ERR : $sample BLAST Failed!" >> $log
                exit 1
            fi
        #output unique contigs if all goes well.
            #flock temporarily locks writing access to a file
            #this makes the parellel tasks take turns writing to the output
                #which stops them from overwriting each other if they finish at the same time
            ( flock -x 9 
            
                echo "sampleID:${sample}" >> ${outdir}/mitotype.out
                #collect all blast ouptut, sort by e-value, report top 10 hits to the output
                cat out/*_blast.out | sort -k 12 -g -r | awk '!a[$12]++' | head -n 10  >> ${outdir}/mitotype.out
            
            ) 9> ${outdir}/.writelock
            
            echo "    $sample FINISHED" >> $log
            #remove temporary files to stop bloat 
            rm -rf tmp/* #lib/*
    done

echo "FINISHED mitotypes: $nstart - $nend" >> $log

date






    

