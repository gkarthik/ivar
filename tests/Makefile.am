LIBS = -lhts -lz -lpthread

CFLAGS = -g -std=c++11

TESTS = check_consensus
check_PROGRAMS = check_consensus
check_consensus_SOURCES = test_call_consensus_from_plup.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
ivar_LDADD = $(LIBS)
