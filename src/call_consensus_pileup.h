#include<iostream>
#include "allele_functions.h"

#ifndef call_consensus_from_pileup
#define call_consensus_from_pileup

void format_alleles(std::vector<allele> &ad);
std::string get_consensus_allele(std::vector<allele> ad, uint8_t min_qual);
int call_consensus_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual);

#endif
