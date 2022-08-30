#include <iostream>
#include <vector>
#include "../src/trim_primer_quality.h"
#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "htslib/sam.h"
#include "../src/primer_bed.h"

int main(){
  int num_tests = 1;
  int success = 0;
  std::vector<primer>::iterator it;
  std::vector<primer> primers = populate_from_file("../data/primer_pair_test/test_primer_pair.bed");
  populate_pair_indices(primers, "../data/primer_pair_test/test_primer_pair_trim.tsv");
  std::string bam = "../data/primer/pair_test/test_primer_pair_trim.untrimmed.bam";    


  //success += flag;
  return (num_tests == success) ? 0 : -1;
}
