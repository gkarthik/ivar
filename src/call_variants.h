#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <regex>
#include <cmath>
#include <htslib/kfunc.h>

#include "allele_functions.h"
#include "ref_seq.h"

#ifndef call_variants
#define call_variants

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double min_threshold, uint8_t min_depth, std::string ref_path, std::string gff_path);
std::vector<allele>::iterator get_ref_allele(std::vector<allele> &ad, char ref);

#endif
