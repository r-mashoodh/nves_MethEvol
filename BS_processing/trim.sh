#!/bin/bash
#!/bin/bash
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -o trim.out
#SBATCH -e trim.err
#SBATCH -J trim


#mkdir software
#git clone https://github.com/FelixKrueger/TrimGalore.git
#cd ..

## trimming bases

mkdir processed_reads
cd raw_data

for i in *_1.fq.gz;
do
filename=`echo $i | awk -F '_1.fq.gz' '{print $1}'`
../software/TrimGalore/trim_galore --paired --fastqc --paired --clip_r1 10 --clip_r2 20 -o ../processed_reads $filename"_1.fq.gz" $filename"_2.fq.gz"
done

# trim_galore --paired --retain_unpaired --illumina --clip_R1 $TrimmedBases --clip_R2 $TrimmedBases --three_prime_clip_R1 $TrimmedBases --three_prime_clip_R2 $TrimmedBases $LeftReads $RightReads


