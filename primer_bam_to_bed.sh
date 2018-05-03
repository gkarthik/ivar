#!/bin/bash

c=$1
p=$2
r=$3

cut -f 1,2 -d ',' $c | sed 's/^/>/g' | tr ',' '\n' > ${p}.fa
~/Documents/bwa/bwa mem -k 5 -T 0 $r  ${p}.fa | samtools view -bS -F 4 -o ${p}.bam
bedtools bamtobed -i ${p}.bam > ${p}.bed
