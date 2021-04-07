#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o bis2.out
#SBATCH -e bis2.err
#SBATCH -J bis2

mkdir bismark_data2
cd bis_aligned

for i in *.deduplicated.bam;
do
../software/Bismark-0.22.3/bismark_methylation_extractor -o ../bismark_data2 --bedGraph --cutoff 10 --zero_based \
--parallel 8 --scaffold --cytosine_report --genome_folder ../genomes/n_vesp/ --gzip $i
done

#mv files to me_data

