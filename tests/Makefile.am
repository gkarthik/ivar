LIBS = -lhts -lz -lpthread

CFLAGS = -g -std=c++11

TESTS = check_consensus check_allele_depth check_consensus_threshold
check_PROGRAMS = check_consensus check_allele_depth check_consensus_threshold
check_consensus_SOURCES = test_call_consensus_from_plup.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
check_allele_depth_SOURCES = test_allele_depth.cpp ../src/allele_functions.cpp
check_consensus_threshold_SOURCES = test_consensus_threshold.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
