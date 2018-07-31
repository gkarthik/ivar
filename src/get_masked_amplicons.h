#include "primer_bed.h"

#ifndef get_masked_amplicons
#define get_masked_amplicons

int get_primer_indice(std::vector<primer> p, unsigned int pos);
int get_primers_with_mismatches(std::string bed, std::string vpath, std::string out);

#endif
