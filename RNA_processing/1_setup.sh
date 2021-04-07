#!/bin/bash

# Setup for new 'Tuxedo' pipeline

mkdir software
cd software

# Download TrimGalore
wget https://github.com/FelixKrueger/TrimGalore/archive/0.5.0.zip
unzip 0.5.0.zip


# Download hisat2
wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip


# Download/install  samtools
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 
tar jxvf samtools-1.9.tar.bz2

cd samtools-1.9
./configure --prefix=`pwd`
make
make install


# Download/install stringtie

git clone https://github.com/gpertea/stringtie

cd stringtie
make release

# Download prepDE.py 
# for creating count matrix with stringtie data
wget http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
chmod +x prepDE.py
