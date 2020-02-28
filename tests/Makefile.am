LIBS = -lhts -lz -lpthread

CXXFLAGS = -g -std=c++11

TESTS = check_quality_trim check_consensus check_allele_depth check_consensus_threshold check_consensus_min_depth check_primer_bed check_getmasked check_removereads check_variants check_common_variants
check_PROGRAMS = check_quality_trim check_consensus check_allele_depth check_consensus_threshold check_consensus_min_depth check_primer_bed check_getmasked check_removereads check_variants check_common_variants
check_quality_trim_SOURCES = check_quality_trim.cpp ../src/trim_primer_quality.cpp ../src/primer_bed.cpp
check_consensus_SOURCES = test_call_consensus_from_plup.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
check_allele_depth_SOURCES = test_allele_depth.cpp ../src/allele_functions.cpp
check_consensus_threshold_SOURCES = test_consensus_threshold.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
check_consensus_min_depth_SOURCES = test_consensus_min_depth.cpp ../src/call_consensus_pileup.cpp ../src/allele_functions.cpp
check_variants_SOURCES = test_variants.cpp ../src/call_variants.cpp ../src/allele_functions.cpp
check_common_variants_SOURCES = test_common_variants.cpp ../src/get_common_variants.cpp
check_primer_bed_SOURCES = test_primer_bed.cpp ../src/primer_bed.cpp
check_getmasked_SOURCES = test_getmasked.cpp ../src/get_masked_amplicons.cpp ../src/primer_bed.cpp
check_removereads_SOURCES = test_removereads.cpp ../src/remove_reads_from_amplicon.cpp ../src/primer_bed.cpp ../src/trim_primer_quality.cpp
