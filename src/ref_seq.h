#include <iostream>
#include <sstream>

#include <htslib/faidx.h>
#include "parse_gff.h"

#ifndef ref_seq
#define ref_seq

class ref_antd{
public:
  ref_antd(std::string ref_path, std::string gff_path);
  char get_base(int64_t pos, std::string region);
  int add_gff(std::string path);
  int add_seq(std::string path);
  std::ostringstream codon_aa_stream(std::ostringstream &line_stream, std::ofstream &fout, int64_t pos, char alt);
  char* get_codon(int64_t pos, std::string region, gff3_feature feature);
  char* get_codon(int64_t pos, std::string region, gff3_feature feature, char alt);

private:
  char *seq;
  gff3 gff;
  faidx_t *fai;
  int ref_len;
  std::string region;
};

#endif
