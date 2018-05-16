build:
	g++ -L/usr/local/lib/ -Lsamtools-1.7/ -lhts main.cpp primer_bed.cpp -lz -lpthread -o trim_primer_quality
	g++ -L/usr/local/lib/ -Lsamtools-1.7/ -lhts remove_reads_from_amplicon.cpp -lz -lpthread -o remove_reads_from_amplicon
