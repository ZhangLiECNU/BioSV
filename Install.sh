#!/bin/bash

#Download tools
#wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
#wget http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2
#wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2
#git clone https://github.com/arq5x/bedtools2.git
#https://codeload.github.com/arq5x/bedtools2/tar.gz/v2.26.0

#Move the downloaded source code into directory:./dependencies/
current_dir=`pwd`

#Unzip tools
tar xfj dependencies/samtools-1.3.1.tar.bz2 -C dependencies/
tar xfj dependencies/bwa-0.7.12.tar.bz2 -C dependencies/
tar zxvf dependencies/bedtools2.tar.gz -C dependencies/
unzip dependencies/wgsim-master.zip -d dependencies/

#Compile Makefile for all tools' makefiles
cd ${current_dir}/dependencies/samtools-1.3.1
./configure --without-curses
make
cd ${current_dir}/dependencies/bwa-0.7.12
make
cd ${current_dir}/dependencies/bedtools2
make
cd ${current_dir}/dependencies/wgsim-master
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm

#copy executable file to tools
cp -p ${current_dir}/dependencies/samtools-1.3.1/samtools ${current_dir}/tools
cp -p ${current_dir}/dependencies/bwa-0.7.12/bwa ${current_dir}/tools
cp -p ${current_dir}/dependencies/bedtools2/bin/bedtools ${current_dir}/tools
cp -p ${current_dir}/dependencies/wgsim-master/wgsim ${current_dir}/tools

cd $current_dir

tools_dir="tools_dir=\""${current_dir}"/tools"\"
Rscript=`which Rscript`
Rscript_header="#!"${Rscript}
Rscript="rscript=\""${Rscript}""\"
echo $Rscript_header > config.txt
echo $tools_dir >> config.txt
echo $Rscript >> config.txt

cat config.txt BioSV.main.R >BioSV
chmod +x BioSV
rm -f config.txt


