#!/bin/bash

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o bis.out
#SBATCH -e bis.err
#SBATCH -J bis

# map bismark

mkdir bis_aligned

cd processed_reads
for i in *_1_val_1.fq.gz; 
do 
filename=`echo $i | awk -F '_1_val_1.fq.gz' '{print $1}'`
../software/Bismark-0.22.3/bismark —bowtie2 --genome ../genomes/n_vesp/ -o ../bis_aligned --score_min L,0,-0.6 -1 $filename"_1_val_1.fq.gz" -2 $filename"_2_val_2.fq.gz"
done

cd ../bis_aligned

for i in *.bam;
do
filename=`echo $i | awk -F '_1_val_1_bismark_bt2_pe.bam' '{print $1}'`
../software/Bismark-0.22.3/deduplicate_bismark $i
done

