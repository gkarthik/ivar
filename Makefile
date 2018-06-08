 # -Lsamtools-1.7/
build:
	g++ -L/usr/local/lib/ -lhts main.cpp primer_bed.cpp -lz -lpthread -o trim_primer_quality
	g++ -L/usr/local/lib/ -lhts remove_reads_from_amplicon.cpp -lz -lpthread -o remove_reads_from_amplicon
	g++ -g -L/usr/local/lib/ -lhts call_consensus.cpp -lz -lpthread -o call_consensus
	g++ -g -L/usr/local/lib/ call_consensus_pileup.cpp -lz -lpthread -o call_consensus_pileup
	g++ -g -L/usr/local/lib/ call_variants.cpp -lz -lpthread -o call_variants
