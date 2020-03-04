#include "../src/ref_seq.h"

int main() {
  int num_success = 0;
  ref_antd refantd("../data/db/test_ref.fa", "../data/test.gff");
  std::vector<gff3_feature> g = refantd.get_gff_features();
  char *codon = new char[3];
  codon = refantd.get_codon(30, "test", g.at(2));
  if(codon[0] != 'C' || codon[1] != 'A' || codon[2] != 'T')
    num_success = -1;
  codon = refantd.get_codon(30, "test", g.at(3));
  if(codon[0] != 'T' || codon[1] != 'C' || codon[2] != 'A')
    num_success = -1;
  codon = refantd.get_codon(29, "test", g.at(2));
  if(codon[0] != 'C' || codon[1] != 'A' || codon[2] != 'T')
    num_success = -1;
  // After edit_pos
  codon = refantd.get_codon(107, "test", g.at(4));
  if(codon[0] != 'C' || codon[1] != 'A' || codon[2] != 'A')
    num_success = -1;
  codon = refantd.get_codon(105, "test", g.at(4));
  if(codon[0] != 'C' || codon[1] != 'A' || codon[2] != 'T')
    num_success = -1;
  codon = refantd.get_codon(101, "test", g.at(4));
  if(codon[0] != 'T' || codon[1] != 'C' || codon[2] != 'T')
    num_success = -1;
  codon = refantd.get_codon(100, "test", g.at(4));
  if(codon[0] != 'G' || codon[1] != 'G' || codon[2] != 'A')
    num_success = -1;
  std::cout << num_success;
  if(num_success == 0)
    return 0;
  return -1;
}
