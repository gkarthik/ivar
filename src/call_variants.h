#include<iostream>

#ifndef call_variants
#define call_variants

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double min_threshold);
std::vector<allele>::iterator get_ref_allele(std::vector<allele> ad, char ref);

#endif
