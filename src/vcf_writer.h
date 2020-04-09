#include <iostream>
#include <fstream>
#include <vector>
#include "allele_functions.h"
#include "ref_seq.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"

class vcf_writer{
  vcfFile *file;
  ref_antd *ref;
  bcf_hdr_t *hdr;
  std::string region;
  std::string sample_name;
  int init_header();
public:
  vcf_writer(char _mode, std::string fname, std::string region, std::string sample_name, std::string ref_path);
  int write_record(uint32_t pos, std::vector<allele> aalt, allele aref, uint32_t depth);
  ~vcf_writer();
};
