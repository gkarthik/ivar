 # -Lsamtools-1.7/
build:
	# g++ -L/usr/local/lib/ -lhts trim_primer_quality.cpp primer_bed.cpp -lz -lpthread -o trim_primer_quality
	# g++ -L/usr/local/lib/ -lhts remove_reads_from_amplicon.cpp -lz -lpthread -o remove_reads_from_amplicon
	# g++ -g -L/usr/local/lib/ -lhts call_consensus.cpp -lz -lpthread -o call_consensus
	# g++ -g -L/usr/local/lib/ call_consensus_pileup.cpp allele_functions.cpp -lz -lpthread -o call_consensus_pileup
	# g++ -g -L/usr/local/lib/ allele_functions.cpp call_variants.cpp -lz -lpthread -o call_variants
	g++ -g alignment.cpp -o alignment
	# g++ -g -L/usr/local/lib/ ivar.cpp call_consensus_pileup.cpp trim_primer_quality.cpp remove_reads_from_amplicon.cpp call_variants.cpp primer_bed.cpp allele_functions.cpp -lhts -lz -lpthread -o ivar
