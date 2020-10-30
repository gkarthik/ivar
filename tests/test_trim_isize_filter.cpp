#include<iostream>
#include "../src/trim_primer_quality.h"
#include "../src/primer_bed.h"
#include "htslib/sam.h"

typedef struct {
  uint16_t flag;
  uint32_t n_cigar;
  hts_pos_t isize;
  uint8_t fail;
} result_t;

int test_trim_isize(uint8_t min_qual, uint8_t sliding_window, bool no_write_flag, bool keep_for_reanalysis, int min_length)
{
  int success = 0;
  int res;

  // bam1_t *aln = bam_init1();

  std::string bam = "../data/test_isize_filter.sorted.bam";
  std::string bed = "../data/test.bed";
  std::string prefix = "trim_isize_filter";
  std::string bam_out = "trim_isize_filter.bam";

  std::string region_ = "";
  std::string cmd = "@PG\tID:ivar-trim\tPN:ivar\tVN:1.0.0\tCL:ivar trim\n";

  // Test and check result
  res = trim_bam_qual_primer(bam, bed, prefix, region_, min_qual, sliding_window, cmd, no_write_flag, keep_for_reanalysis, min_length);
  return success;
}
int main() {
  int success = 0;

  // Default paramaters
  success = test_trim_isize(20, 4, true, false, 30);

  return success;
}


