#!/bin/bash

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
