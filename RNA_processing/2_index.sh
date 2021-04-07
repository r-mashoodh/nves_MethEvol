#!/bin/bash

# Build Hisat2 index

# Download beetle genome fasta
mkdir genome
cd genome

wget  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/412/225/GCF_001412225.1_Nicve_v1.0/GCF_001412225.1_Nicve_v1.0_genomic.fna.gz
gunzip GCF_001412225.1_Nicve_v1.0_genomic.fna.gz


# Download beetle gtf
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Nicrophorus_vespilloides/GFF/ref_Nicve_v1.0_top_level.gff3.gz
gunzip ref_Nicve_v1.0_top_level.gff3.gz

# Extract splice-site and exon information from the gene annotation file
cd ..
mkdir hisat2_index
cd hisat2_index

../software/hisat2-2.1.0/extract_splice_sites.py ../genome/ref_Nicve_v1.0_top_level.gff3 >beetle.ss
../software/hisat2-2.1.0/extract_exons.py ../genome/ref_Nicve_v1.0_top_level.gff3 >beetle.exon


# Build index
../software/hisat2-2.1.0/hisat2-build --ss beetle.ss --exon beetle.exon ../genome/GCF_001412225.1_Nicve_v1.0_genomic.fna beetle_index
