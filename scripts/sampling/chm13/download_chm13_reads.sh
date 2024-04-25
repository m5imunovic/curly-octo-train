#!/bin/bash
READS_DIR="./fastq"
mkdir $READS_DIR

wget -O SRX7897685.fastq https://sra-pub-src-2.s3.amazonaws.com/SRR11292123/m64062_190804_172951.fastq.1
mv SRX7897685.fastq $READS_DIR
wget -O SRX7897686.fastq https://sra-pub-src-2.s3.amazonaws.com/SRR11292122/m64062_190807_194840.fastq.1
mv SRX7897686.fastq $READS_DIR
wget -O SRX7897687.fastq https://sra-pub-src-2.s3.amazonaws.com/SRR11292121/m64062_190803_042216.fastq.1
mv SRX7897687.fastq $READS_DIR
wget -O SRX7897688.fastq https://sra-pub-src-2.s3.amazonaws.com/SRR11292120/m64062_190806_063919.fastq.1
mv SRX7897688.fastq $READS_DIR
