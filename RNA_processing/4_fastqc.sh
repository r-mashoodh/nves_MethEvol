#!/bin/bash

# fastqc for adaptor trimmed reads

mkdir fastqc_files
cd processed_reads

for i in *.fq.gz;
do
filename =`echo $i | awk -F '.fq.gz' '{print $1}'`
fastqc $filename".fq.gz" -o ../fastqc_files
done

