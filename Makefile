build:
	g++ -L/usr/local/lib/ -Lsamtools-1.7/ -lhts main.cpp -lz -lpthread -o ivar
