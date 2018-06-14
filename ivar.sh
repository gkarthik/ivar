#!/bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
usage_dir=$script_dir/usage

usage() {
    exec 1>&2
    echo
    case "$1" in
	trim)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	removereads)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	callvariants)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	filtervariants)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	consensus)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	createbed)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	getmasked)
	    cat ${usage_dir}/usage_$1.txt;
	    ;;
	*)
	    cat ${usage_dir}/usage_all.txt;
	    ;;
    esac
    echo
    exit 1;
}

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
		    usage trim
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ] || [ -z "${b}" ]; then
	    usage trim
	fi
	$script_dir/trim_primer_quality ${i} ${b} ${p}.trimmed.bam
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
		    usage removereads
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ]; then
	    usage removereads
	fi
	$script_dir/remove_reads_from_amplicon ${i} ${p}.filtered.bam "$@"
	echo "Sorting"
	samtools sort -o ${p}.filtered.sorted.bam ${p}.filtered.bam
	echo "Indexing"
	samtools index ${p}.filtered.sorted.bam ${p}.filtered.sorted.bai
	;;
    callvariants)
	while getopts ":i:p:r:R:q:" o; do
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
		R)
		    R=${OPTARG}
		    ;;
		q)
		    q=${OPTARG}
		    ;;
		*)
		    usage callvariants
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ] || [ -z "${r}" ]; then
	    usage callvariants
	fi
	if [ -z "$R" ]; then
	    samtools mpileup -A -B -Q 0 -d 300000 -pm 1 -F 0 --reference ${r} ${i} | $script_dir/call_variants ${p} ${q:-20}
	else
	    samtools mpileup -A -B -Q 0 -d 300000 -pm 1 -F 0 -r ${R} --reference ${r} ${i} | $script_dir/call_variants ${p} ${q:-20}
	fi
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
		*)
		    usage filtervariants
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${p}" ] || [ -z "${f}" ]; then
	    usage filtervariants
	fi
	$script_dir/combine_variants.py ${p} ${f} "$@"
	;;
    getmasked)
	while getopts ":p:f:" o; do
	    case "${o}" in
		p)
		    p=${OPTARG}
		    ;;
		f)
		    f=${OPTARG}
		    ;;
		*)
		    usage getmasked
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${p}" ] || [ -z "${f}" ]; then
	    usage getmasked
	fi
	$script_dir/get_masked_amplicons.py ${p} ${f} "$@"
	;;
    consensus)
	while getopts ":i:p:r:R:q:" o; do
	    case "${o}" in
		i)
		    i=${OPTARG}
		    ;;
		r)
		    r=${OPTARG}
		    ;;
		p)
		    p=${OPTARG}
		    ;;
		R)
		    R=${OPTARG}
		    ;;
		q)
		    q=${OPTARG}
		    ;;
		*)
		    usage consensus
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${i}" ] || [ -z "${p}" ]; then
	    usage consensus
	fi
	if [ -z "$R" ]; then
	    samtools mpileup -A -B -Q 0 -d 300000 -pm 1 -F 0 --reference ${r} ${i} | $script_dir/call_consensus_pileup ${p} ${q:-20}
	else
	    samtools mpileup -A -B -Q 0 -d 300000 -pm 1 -F 0 -r ${R} --reference ${r} ${i} | $script_dir/call_consensus_pileup ${p} ${q:-20}
	fi
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
		    usage createbed
		    ;;
	    esac
	done
	shift $((OPTIND-1))
	if [ -z "${c}" ] || [ -z "${r}" ] || [ -z "${p}" ]; then
	    usage createbed
	fi
	$script_dir/primer_bam_to_bed.sh $c $p $r
	;;
    *)
	usage
esac
