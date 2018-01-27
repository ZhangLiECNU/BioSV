#!/bin/bash
if [ $# != 1 ]
then
        echo "Usege:sh $0 <UCSC reference fasta file>"
        echo "please input the Reference genome file, from which only chromosome chr21 and chr22 will be chosen for further analysis."
        exit
fi

Ref=$1
#echo "before run this example please chech if the chromosomes names in chromosomes.txt are in accordance with the chromosomes names in your reference file (please note 'chr' in chromosomes.txt)"

#Generating a fasta file that only include chr21 and chr22
#./tools/samtools faidx ${Ref} chr21 chr22 > example.fa
cat ./test_data/chromosomes.txt |xargs ./tools/samtools faidx ${Ref} > test_data/test.fa

#Build bwa index for example_N.fa, which only include chromosomes of chr21 and chr22
./tools/bwa index test_data/test.fa

#Simulating structural variation from diploid genomes (chr21 and chr22)
python2.7 BioSV_simulator.py -g ./test_data/test.fa -o test_data/simulate -p example -c ./test_data/chromosomes.txt -n 50,0,50,50,50 -N 10000000

### data preprocessing
./tools/bwa mem ./test_data/test.fa ./test_data/simulate/example_N1.fq ./test_data/simulate/example_N2.fq | ./tools/samtools view -bS - > ./test_data/test.bam
./tools/samtools sort -o ./test_data/test.sort.bam ./test_data/test.bam
./tools/samtools index ./test_data/test.sort.bam
#Structural variation calling
./BioSV call -b ./test_data/test.sort.bam -o test_data/BioSV_output -t 2 --minMRC 4 

awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,NR,".",".",".",$9,$10}' ./test_data/simulate/example_N.bedpe >./test_data/simulate/simu.bedpe
grep -v "#" ./test_data/BioSV_output/BioSV.output.bedpe | grep  "PASS" | awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,NR,".",".",".",$7,substr($9,4,3)}' > ./test_data/BioSV_output/predicted.hc.bedpe
Rscript ./test_data/performance.R

