#!/bin/bash

#mkdir genomes
#cd genomes

#mkdir n_vesp
#wget  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/412/225/GCF_001412225.1_Nicve_v1.0/GCF_001412225.1_Nicve_v1.0_genomic.fna.gz
#gunzip GCF_001412225.1_Nicve_v1.0_genomic.fna.gz

cd genomes 
# genome preparation for m_zebra
../software/Bismark-0.22.3/bismark_genome_preparation --bowtie2 n_vesp/

# genome preparation for spike in
#../software/Bismark-0.22.3/bismark_genome_preparation --path_to_bowtie ../software/bowtie2 spike/
