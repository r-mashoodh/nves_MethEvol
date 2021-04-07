#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o requant.out
#SBATCH -e requant.err
#SBATCH -J requant

# execute within requant dir
# this script takes aligned files from hisat and uses stringtie to quantify genes using NCBIs annotation
# rather than a faux de novo based on novel transcripts identified by hisat

mkdir ballgown

cd gtf_files

for i in *.gtf
do
filename=`echo $i | awk -F ".gtf" '{print $1}'`
../../software/stringtie/stringtie -e -B -p 8 -G ../../genome/ref_Nicve_v1.0_top_level.gff3 -o ../ballgown/$filename/$filename".gtf" ../../sam_files/$filename".sorted.bam"
done

# create table using prepDE.py (from stringtie)

cd ..

../software/prepDE.py
