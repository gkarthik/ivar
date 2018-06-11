#!/bin/bash

# ./call_variants.sh temp.sorted.bam ~/hpc_downloads/2018.02.27/FSS13025.fasta calls.vcf.gz temp.consensus.fa

in=$1
ref=$2
out=$3
samtools mpileup -m 1 -A -B -Q 0 -d 300000 -Ou -t AD -f $ref $in | bcftools call -m -A -Oz -o ./variants/$out.vcf.gz
bcftools index $out.vcf.gz

# date && find ./ -name "*.sorted.bam" -exec basename {} \; | sed 's/.bam//g' | xargs -I {} -L 1 -P 8 sh -c "samtools mpileup -A -B -Q 0 -d 300000 -Ou -t AD -f ~/Documents/code/ivar/data/db/ZIKV_PRV.fasta {}.bam | bcftools call --ploidy 1 -m -A -Oz -o ./variants/{}.vcf.gz" && date
