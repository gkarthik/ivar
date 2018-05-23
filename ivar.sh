#!/bin/bash

usage() { echo "Usage: $0 [command <trim|callvariants|filtervariants|consensus|createbed|removereads|getmasked>]" 1>&2; exit 1; }
usage_trim() { echo "Usage: $0 trim [-i input] [-b <bedfile>] [-p prefix]" 1>&2; exit 1; }
usage_removereads() { echo "Usage: $0 removereads [-i input] [-p prefix] <amplicon-indice-1> <amplicon-indice-2> ..." 1>&2; exit 1; }
usage_call_variants() { echo "Usage: $0 callvariants [-i input] [-r reference] [-p prefix]" 1>&2; exit 1; }
usage_filter_variants() { echo "Usage: $0 filtervariants [-f <frequency cut off>] [-b <bed file>] [-p prefix] replicate1.vcf.gz replicate2.vcf.gz ... " 1>&2; exit 1; }
usage_consensus() { echo "Usage: $0 consensus [-i <input-vcf>] [-p prefix] [-r reference] " 1>&2; exit 1; }
usage_create_bed() { echo "Usage: $0 createbed [-c <primer-csv>] [-p prefix] [-r reference] " 1>&2; exit 1; }
usage_get_masked() { echo "Usage: $0 getmasked [-f <frequency-cut-off>] [-p prefix] [-b bed-file] replicate1.vcf.gz replicate2.vcf.gz ..." 1>&2; exit 1; }

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
    removereads)
	while getopts ":i:p:a:" o; do
	    case "${o}" in
		i)
		    i=${OPTARG}
		    ;;
		p)
		    p=${OPTARG}
		    ;;
		*)
		    usage_removereads
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ]; then
	    usage_removereads
	fi
	~/Documents/code/ivar/remove_reads_from_amplicon ${i} ${p}.filtered.bam "$@"
	echo "Sorting"
	samtools sort -o ${p}.filtered.sorted.bam ${p}.filtered.bam
	echo "Indexing"
	samtools index ${p}.filtered.sorted.bam ${p}.filtered.sorted.bai
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
	samtools mpileup -A -B -Q 0 -d 300000 -gu -t AD -pm 1 -F 0 -f ${r} ${i} | bcftools call --ploidy 1 -m -A -Oz -o ${p}.vcf.gz
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
    getmasked)
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
		    usage_get_masked
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${p}" ] || [ -z "${f}" ] || [ -z "${b}" ]; then
	    usage_get_masked
	fi
	~/Documents/code/ivar/get_masked_amplicons.py ${p} ${f} ${b} "$@"
	;;
    consensus)
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
		    usage_consensus
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${r}" ] || [ -z "${p}" ]; then
	    usage_consensus
	fi
	
	cat ${r} | bcftools consensus $i > ${p}.fa
	;;
    createbed)
	while getopts ":c:p:r:" o; do
	    case "${o}" in
		c)
		    c=${OPTARG}
		    ;;
		p)
		    p=${OPTARG}
		    ;;
		r)
		    r=${OPTARG}
		    ;;
		*)
		    usage_create_bed
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${r}" ] || [ -z "${p}" ]; then
	    usage_create_bed
	fi
	~/Documents/code/ivar/primer_bam_to_bed.sh $c $p $r
	;;
    *)
	usage
esac
