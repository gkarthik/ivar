#!/bin/bash

usage() { echo "Usage: $0 [command <trim|callvariants|filtervariants|consensus>]" 1>&2; exit 1; }
usage_trim() { echo "Usage: $0 trim [-i input] [-b <bedfile>] [-p prefix]" 1>&2; exit 1; }
usage_call_variants() { echo "Usage: $0 callvariants [-i input] [-r reference] [-p prefix]" 1>&2; exit 1; }
usage_filter_variants() { echo "Usage: $0 filtervariants [-f <frequency cut off>] [-b <bed file>] [-p prefix] replicate1.vcf.gz replicate2.vcf.gz ... " 1>&2; exit 1; }

cmd=$1; shift
case "$cmd" in
    trim)
	while getopts ":i:p:b:" o; do
	    case "${o}" in
		i)
		    i=${OPTARG}
		    ;;
		b)
		    b=${OPTARG}
		    ;;
		p)
		    p=${OPTARG}
		    ;;
		*)
		    usage_trim
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ] || [ -z "${b}" ]; then
	    usage_trim
	fi
	~/Documents/code/ivar/trim_primer_quality ${i} ${b} ${p}.trimmed.bam
	echo "Sorting"
	samtools sort -o ${p}.trimmed.sorted.bam ${p}.trimmed.bam
	echo "Indexing"
	samtools index ${p}.trimmed.sorted.bam ${p}.trimmed.sorted.bai
	;;
    callvariants)
	while getopts ":i:p:r:" o; do
	    case "${o}" in
		i)
		    i=${OPTARG}
		    ;;
		p)
		    p=${OPTARG}
		    ;;
		r)
		    r=${OPTARG}
		    ;;
		*)
		    usage_call_variants
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ] || [ -z "${r}" ]; then
	    usage_call_variants
	fi
	samtools mpileup -A -B -Q 0 -d 300000 -Ou -t AD -f ${r} ${i} | bcftools call --ploidy 1 -m -A -Oz -o ${p}.vcf.gz
	bcftools index ${p}.vcf.gz
	;;
    filtervariants)
	while getopts ":p:f:b:" o; do
	    case "${o}" in
		p)
		    p=${OPTARG}
		    ;;
		f)
		    f=${OPTARG}
		    ;;
		b)
		    b=${OPTARG}
		    ;;
		*)
		    usage_filter_variants
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${p}" ] || [ -z "${f}" ] || [ -z "${b}" ]; then
	    usage_filter_variants
	fi
	~/Documents/code/ivar/combine_variants.py ${p} ${f} ${b} "$@"
	;;
    consensus)
	while getopts ":i:p:" o; do
	    case "${o}" in
		i)
		    i=${OPTARG}
		    ;;
		p)
		    p=${OPTARG}
		    ;;
		*)
		    usage
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ]; then
	    usage
	fi
	;;
    *)
	usage
	;;
esac
