#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o 3_processing.out
#SBATCH -e 3_processing.err
#SBATCH -J 3_processing

# filenames eg.
#N1_76_1.fq
#N1_76_2.fq
#F1_09a_1.fq
#F1_09a_2.fq

## ADAPTOR TRIMMING ##

mkdir processed_reads
cd raw_data

for i in *_1.fq.gz;
do
filename=`echo $i | awk -F '_1.fq.gz' '{print $1}'`
../software/TrimGalore-0.5.0/trim_galore --illumina --paired -o ../processed_reads $filename"_1.fq.gz" $filename"_2.fq.gz"
done

cd ../processed_reads

## FASTQC ##

# ???

#

## ALIGNMENT WITH HISAT2

# files are contained in processed_reads dir 
# file names should/will be sample_1_val_1.fq etc

mkdir ../sam_files

for i in *_1_val_1.fq.gz; 
do 
filename=`echo $i | awk -F '_1_val_1.fq.gz' '{print $1}'`
../software/hisat2-2.1.0/hisat2 -p 8 --dta -x ../hisat2_index/beetle_index -1 $filename"_1_val_1.fq.gz" -2 $filename"_2_val_2.fq.gz" -S ../sam_files/$filename.sam
done


# Sort and convert the SAM  les to BAM note: command works with SAMtools version 1.3 or newer

cd ../sam_files

for i in *.sam
do
filename=`echo $i | awk -F "." '{print $1}'`
../software/samtools-1.9/samtools view -bS $filename.sam -o $filename.bam
../software/samtools-1.9/samtools sort $filename.sam -o $filename.sorted.bam
done


## ASSEMBLE WITH STRINGTIE
mkdir ../gtf_files

for i in *.sorted.bam
do
filename=`echo $i | awk -F ".sorted.bam" '{print $1}'`
../software/stringtie/stringtie -p 8 -G ../genome/ref_Nicve_v1.0_top_level.gff3 -o ../gtf_files/$filename.gtf -l $filename $filename.sorted.bam
done

cd ../gtf_files

files=("$(ls)")
printf '%s\n' "${files[@]}" > mergelist.txt

#for i in ../gtf_files/*.gtf
#do
#echo $i >> mergelist.txt
#done

# merge all transcripts
../software/stringtie/stringtie --merge -p 8 -G ../genome/ref_Nicve_v1.0_top_level.gff3 -o ../stringtie_merged.gtf mergelist.txt

mkdir ../ballgown

# create table counts assemble with merged gtf using stringtie

for i in *.gtf
do
filename=`echo $i | awk -F ".gtf" '{print $1}'`
../software/stringtie/stringtie -e -B -p 8 -G ../stringtie_merged.gtf -o ../ballgown/$filename/$filename".gtf" ../sam_files/$filename".sorted.bam"
done

# create table using prepDE.py (from stringtie)
cd ..
./software/prepDE.py

# https://www.biostars.org/p/282817/ about the MSTRG geneID

