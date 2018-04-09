#!/bin/bash

# ./call_consensus.sh temp.sorted.bam ~/hpc_downloads/2018.02.27/FSS13025.fasta calls.vcf.gz temp.consensus.fa

in=$1
ref=$2
out=$3
consensus_out=$4
bcftools mpileup -d 1000000 -Ou -f $ref $in | bcftools call --ploidy 1 -mv -Oz -o $out

bcftools index $out
# Consensus calling
samtools faidx $ref KU955593.1 | bcftools consensus $out > $consensus_out
