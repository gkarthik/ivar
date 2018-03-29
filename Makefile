build:
	g++ -L/usr/local/lib/ -Lsamtools-1.7/ -lhts main.cpp primer_bed.cpp -lz -lpthread -o ivar
